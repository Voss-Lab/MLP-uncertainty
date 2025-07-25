"""
ASE's mixing.py calculator is customized by Johannes Voss and Suman Bhasker Ranganath:
    Added multiprocessing for parallel execution of ASE based calculators (ML potentials)
    Used shared memory for efficient storage/access of atomic positions, energy, and forces
    Individual ML ensemble predictions are used to quantify uncertainty and compute average
    Implemented metadynamics based exploration of atomic configurations with forcemode=2,3    
"""

from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
import multiprocessing
from multiprocessing import shared_memory
import numpy as np

def __multiproccalc__(cpipe, calctype, atoms, spositionsname, senergyname, sforcesname, calcargs, kwargs):
    calc = calctype(*calcargs, **kwargs)
    atoms.set_calculator(calc)
    shmpos = shared_memory.SharedMemory(name=spositionsname)
    shmforces = shared_memory.SharedMemory(name=sforcesname)
    shmenergy = shared_memory.SharedMemory(name=senergyname)

    pos = np.ndarray(atoms.positions.shape, dtype=atoms.positions.dtype, buffer=shmpos.buf)
    forces = np.ndarray(atoms.positions.shape, dtype=atoms.positions.dtype, buffer=shmforces.buf)
    energy = np.ndarray((1,), dtype=atoms.positions.dtype, buffer=shmenergy.buf)

    while(cpipe.recv()):
        atoms.set_positions(pos)
        forces[:][:] = atoms.get_forces()
        energy[0] = atoms.get_potential_energy()
        cpipe.send(True)

    del calc
    shmpos.close()
    shmforces.close()
    shmenergy.close()
    cpipe.close()


class MultiProcessCalculator(Calculator):
    def __init__(self, *args, **kwargs):
        if 'atoms' not in kwargs.keys():
            raise ValueError('atoms=AtomsObject must be specified as parameter.')
        calctype = args[0]
        calcargs = args[1:]
        super().__init__(**kwargs)
        calckwargs = kwargs.copy()
        calckwargs.pop('atoms')

        self.atoms = kwargs['atoms']
        self.spositions = shared_memory.SharedMemory(create=True, size=self.atoms.positions.nbytes)
        self.sforces = shared_memory.SharedMemory(create=True, size=self.atoms.positions.nbytes)
        self.senergy = shared_memory.SharedMemory(create=True, size=self.atoms.positions[0][:1].nbytes)

        self.spos = np.ndarray(self.atoms.positions.shape, dtype=self.atoms.positions.dtype, buffer=self.spositions.buf)
        self.sfor = np.ndarray(self.atoms.positions.shape, dtype=self.atoms.positions.dtype, buffer=self.sforces.buf)
        self.sene = np.ndarray((1,), dtype=self.atoms.positions.dtype, buffer=self.senergy.buf)

        self.implemented_properties = {
            'energy' : self.calculate_energy_and_forces,
            'forces' : self.calculate_energy_and_forces,
        }
        self.results = {}

        self.ppipe, self.cpipe = multiprocessing.Pipe()
        self.process = multiprocessing.Process(target=__multiproccalc__, args=(self.cpipe,calctype,self.atoms,self.spositions.name,self.senergy.name,self.sforces.name,calcargs,calckwargs))
        self.process.start()

    def stop(self):
        self.ppipe.send(False)
        self.process.join()
        self.spositions.close()
        self.sforces.close()
        self.senergy.close()

    def __del__(self):
        self.stop()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        super().calculate(atoms, properties, system_changes)

        has_results = [p in self.results for p in properties]
        if (len(system_changes) > 0) or (not np.all(has_results)):
            if len(atoms)!=len(self.atoms) or np.product(atoms.numbers==self.atoms.numbers)==0 or np.product(atoms.cell==self.atoms.cell)==0:
                raise ValueError('MultiProcessCalculator cannot deal with atomic species or unit cell changes.')
            self.spos[:][:] = atoms.positions
            if ('energy' in properties) and ('forces' in properties):
                del properties[properties.index('energy')]
            for p in properties:
                if p in self.implemented_properties:
                    self.implemented_properties[p](self.atoms)
                else:
                    raise NotImplementedError(
                        "Property not implemented: {}".format(p))

    #def calculate_energy(self, atoms):
    #    self.ppipe.send(True)
    #    self.results['energy'] = 0.0

    def calculate_energy_and_forces(self, atoms):
        self.ppipe.send(True)
        self.results['energy'] = 0.0
        self.results['forces'] = self.sfor

    def multiprocsync(self):
        if not self.ppipe.recv():
            raise RuntimeError('MultiProcessCalculator subprocess failed.')
        self.results['energy'] = self.sene[0]


class LinearCombinationCalculator(Calculator):
    """LinearCombinationCalculator for weighted summation of multiple calculators.
    """

    def __init__(self, calcs, weights, atoms=None, forcemode='average', varscale=0.05, metascale=0.5):
        """Implementation of sum of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        weights: list of float
            Weights for each calculator in the list.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        super().__init__(atoms=atoms)

        if len(calcs) == 0:
            raise ValueError('The value of the calcs must be a list of Calculators')

        for calc in calcs:
            if not isinstance(calc, Calculator):
                raise ValueError('All the calculators should be inherited form the ase\'s Calculator class')

        common_properties = set.intersection(*(set(calc.implemented_properties) for calc in calcs))
        self.implemented_properties = list(common_properties)

        if not self.implemented_properties:
            raise PropertyNotImplementedError('There are no common property implemented for the potentials!')

        if len(weights) != len(calcs):
            raise ValueError('The length of the weights must be the same as the number of calculators!')

        self.calcs = calcs
        self.weights = weights
        self.forcemode = ['average',
                          'vargrad',
                          'meta',
                          'metamix'].index(forcemode)
        self.varscale = varscale
        self.metascale = metascale

        self.std_energy = None
        self.std_forces = None

    def calculate(self, atoms=None, properties=['energy','forces'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)

        if not set(properties).issubset(self.implemented_properties):
            raise PropertyNotImplementedError('Some of the requested property is not in the '
                                              'list of supported properties ({})'.format(self.implemented_properties))

        waitforcalcs = []

        for w, calc in zip(self.weights, self.calcs):
            if calc.calculation_required(atoms, properties):
                waitforcalcs.append(calc)
                calc.calculate(atoms, properties, system_changes)

        for calc in waitforcalcs:
            calc.multiprocsync()

        properties = ['energy', 'forces']

        for w, calc in zip(self.weights, self.calcs):
            for k in properties:
                if k not in self.results:
                    self.results[k] = w * calc.results[k]
                else:
                    self.results[k] += w * calc.results[k]

        if len(waitforcalcs)>0:
            self.std_energy = 0.0
            self.std_forces = 0.0
            energy = self.results['energy']
            forces = self.results['forces']
            for w, calc in zip(self.weights, self.calcs):
                self.std_energy += w*(calc.results['energy']-energy)**2
                self.std_forces += w*np.linalg.norm(calc.results['forces']-forces)**2
            self.std_energy = self.std_energy**0.5
            self.std_forces = self.std_forces**0.5

            if self.forcemode>0:
                self.Eave = self.results['energy']
                self.Fave = self.results['forces'].copy()
                forces = np.zeros_like(self.Fave)
                for w, calc in zip(self.weights, self.calcs):
                    forces += w*(calc.results['energy']-self.Eave) * (calc.results['forces']-self.Fave)

                forces *= 2./len(forces)
                varnorm = self.std_energy**2/len(forces)

                if self.forcemode==1:
                    self.results['energy'] = varnorm
                    self.results['forces'] = forces
                elif self.forcemode==2:
                   exp = self.metascale * np.exp(-varnorm / self.varscale)
                   self.results['energy'] = self.metascale * np.exp(-varnorm / self.varscale) * len(forces) #bias Energy                 
                   self.results['forces'] = -self.metascale * np.exp(-varnorm / self.varscale) * forces / self.varscale * len(forces)
                elif self.forcemode==3:
                    exp = self.metascale * np.exp(-varnorm / self.varscale)
                    self.results['energy'] = self.Eave + len(forces) * exp
                    self.results['forces'] = self.Fave - len(forces) * exp * forces / self.varscale

    def reset_varscale(self, varscale):
        if self.forcemode<2:
           raise RuntimeError('Variance-resetting only makes sense for metadynamics')
        self.varscale = varscale

    def reset(self):
        """Clear all previous results recursively from all of the calculators."""
        super().reset()

        for calc in self.calcs:
            calc.reset()

    def get_std_energy(self):
        """Returns standard deviation of total energy from ensemble.
           Does not check if energies need to be recomputed."""
        return self.std_energy

    def get_std_forces(self):
        """Returns standard deviation of forces from ensemble.
           Does not check if forces need to be recomputed."""
        return self.std_forces
    
    def get_ensemble_avg_energy(self):
        """Returns average of total energy from ensemble.
           Does not check if energies need to be recomputed."""
        return self.Eave
        
    def __str__(self):
        calculators = ', '.join(calc.__class__.__name__ for calc in self.calcs)
        return '{}({})'.format(self.__class__.__name__, calculators)


class MixedCalculator(LinearCombinationCalculator):
    """
    Mixing of two calculators with different weights

    H = weight1 * H1 + weight2 * H2

    Has functionality to get the energy contributions from each calculator

    Parameters
    ----------
    calc1 : ASE-calculator
    calc2 : ASE-calculator
    weight1 : float
        weight for calculator 1
    weight2 : float
        weight for calculator 2
    """

    def __init__(self, calc1, calc2, weight1, weight2):
        super().__init__([calc1, calc2], [weight1, weight2])

    def set_weights(self, w1, w2):
        self.weights[0] = w1
        self.weights[1] = w2

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """ Calculates all the specific property for each calculator and returns with the summed value.
        """

        super().calculate(atoms, properties, system_changes)
        if 'energy' in properties:
            energy1 = self.calcs[0].get_property('energy', atoms)
            energy2 = self.calcs[1].get_property('energy', atoms)
            self.results['energy_contributions'] = (energy1, energy2)

    def get_energy_contributions(self, atoms=None):
        """ Return the potential energy from calc1 and calc2 respectively """
        self.calculate(properties=['energy'], atoms=atoms)
        return self.results['energy_contributions']


class SumCalculator(LinearCombinationCalculator):
    """SumCalculator for combining multiple calculators.

    This calculator can be used when there are different calculators for the different chemical environment or
    for example during delta leaning. It works with a list of arbitrary calculators and evaluates them in sequence
    when it is required.
    The supported properties are the intersection of the implemented properties in each calculator.
    """

    def __init__(self, calcs, atoms=None):
        """Implementation of sum of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        weights = [1.] * len(calcs)
        super().__init__(calcs, weights, atoms)


class AverageCalculator(LinearCombinationCalculator):
    """AverageCalculator for equal summation of multiple calculators (for thermodynamic purposes)..
    """

    def __init__(self, calcs, atoms=None, **kwargs):
        #print(*kwargs)
        """Implementation of average of calculators.

        calcs: list
            List of an arbitrary number of :mod:`ase.calculators` objects.
        atoms: Atoms object
            Optional :class:`~ase.Atoms` object to which the calculator will be attached.
        """

        n = len(calcs)

        if n == 0:
            raise ValueError('The value of the calcs must be a list of Calculators')

        weights = [1 / n] * n
        super().__init__(calcs, weights, atoms, **kwargs)
