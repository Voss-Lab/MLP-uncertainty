"""
This script analyzes the evolution of local atomic environment 
in MD trajectory of Pt(111)-hydrogen interface, generating 4-panel plots.
Created by Filippo Balzaretti, modified by Johannes Voss and Suman Bhasker-Ranganath
"""

from ase import neighborlist
from ase.io.trajectory import Trajectory
from ase.geometry import get_distances
import numpy as np
import matplotlib.pyplot as plt
import os
from ase.io.jsonio import write_json
import sys



# --------- DICTIONARIES  ---------- #
font = {'family': 'DejaVu Sans',
        'style': 'normal',
        'color':  'black',
        'weight': 'normal',
        'size': 35,
        }

ads_sites = {'0': 'vacuum' ,
             '1': 'ontop'  ,
             '2': 'bridge' ,
             '3': 'hollow' }

site_colors = {'vacuum' : '#1f77b4', # default matplotlib color
               'ontop'  : 'green'  ,
               'bridge' : 'orange' ,
               'hollow' : 'purple' ,
                           'in bulk': 'black'}

itr_colors = {'gas H'     : 'blue',
                          'gas H\u2099'  : 'red'  ,
                          'H*'    : 'teal' ,
                          'H\u2099*' : 'darkred' }



# ---------- FUNCTIONS ----------- #
def parse_input():
        if len(sys.argv) == 2 and sys.argv[1] == '--h':
                print("""
EXAMPLE OF USAGE:
       - python3 stats.py run.traj_all H_all Pt_all

INSTRUCTIONS:
The script runs with user inputs in the following order:
Usage: python3 stats.py <trajectory> <adsorbant atoms> <slab atoms>

<trajectory>:
   Name of trajectory file followed by the frames of interest.
   Use 'all' for considering all frames.
   Examples of indexing for a trajectory file 'run.traj':
       - run.traj_[0:100]
       - run.traj_[50,100,150]
       - run.traj_all

<adsorbant atoms>:
   Name of adsorbant species followed by their indeces of interest.
   Use 'all' for considering all adsorbant indeces.
   Examples on indexing for the adsorbant species 'H':
       - H_[0:10]
       - H_[0,10,19]
       - H_all

<slab atoms>:
   Name of slab species followed by their indeces of interest.
   Use 'all' for considering all slab indeces.
   Examples on indexing for the slab species 'Pt':
       - Pt_[0:10]
       - Pt_[0,10,19]
       - Pt_all

ADDITIONAL INFORMATION:
- Make sure to separate with the underscore '_' and not to have it anywhere else
- Make sure to only have empty spaces ' ' between the input and not anywhere else.
  For example, H_[1,2] works, but H_[1, 2] fails.
""")
                sys.exit()

        elif len(sys.argv) == 4:
                traj      = Trajectory(sys.argv[1].split('_')[0])
                traj_inds = sys.argv[1].split('_')[1]
                ads       = sys.argv[2].split('_')[0]
                ads_inds  = sys.argv[2].split('_')[1]
                slb       = sys.argv[3].split('_')[0]
                slb_inds  = sys.argv[3].split('_')[1]

                if traj_inds == 'all':
                        frames = traj
                else:
                        if ':' not in traj_inds:
                                traj_indeces = eval(traj_inds)
                        else:
                                start, end  = map(int, traj_inds[1:-1].split(':'))
                                traj_indeces = list(range(start, end))

                        frames = [traj[i] for i in traj_indeces]

                if ads_inds == 'all':
                        ads_indeces = [i for i, atom in enumerate(frames[-1])
                                           if atom.symbol == ads]
                else:
                        if ':' not in ads_inds:
                                ads_indeces = eval(ads_inds)

                        else:
                                start, end  = map(int, ads_inds[1:-1].split(':'))
                                ads_indeces = list(range(start, end))

                if slb_inds == 'all':
                        slb_indeces = [i for i, atom in enumerate(frames[-1])
                                           if atom.symbol == slb]
                else:
                        if ':' not in slb_inds:
                                slb_indeces = eval(slb_inds)

                        else:
                                start, end  = map(int, slb_inds[1:-1].split(':'))
                                slb_indeces = list(range(start, end))


                return traj, frames, ads, ads_indeces, slb, slb_indeces
        else:
                print("You need 3 input string. Use '--h' for more details.")
                sys.exit()


def get_connectivity_matrix(atoms_object):
    """
    Generate a connectivity matrix for the atoms_object.

    Parameters:
    - atoms_object: ASE Atoms object
        The atoms object for which the connectivity matrix is to be generated.

    Returns:
    - connectivity_matrix: 2D numpy array
        The connectivity matrix where each element (i, j) represents whether
        atom i is bonded to atom j. A value of 1 indicates a bond, and 0 indicates
        no bond.
    """

    cutoff = neighborlist.natural_cutoffs(atoms_object) # covalent radii cutoff
    neighborList = neighborlist.NeighborList(cutoff, self_interaction=False,
                                                                                         bothways=True)
    neighborList.update(atoms_object)
    connectivity_matrix = neighborList.get_connectivity_matrix(sparse=False)

    return connectivity_matrix

def plot_stats(information):
        """
    Plot statistical information for a single adsorbant in the dynamics.

    Parameters:
    - information: list of dictionaries
        A list where each dictionary represents statistical information for
        a specific frame. Each dictionary should contain keys 'slb_dist',
        'site_type', 'tot_coord', and 'ads_type', corresponding to the
        surface distance, adsorption site type, total coordination, and
        adsorption interaction type, respectively.

    Returns:
    None, but it saves the figure in the Figs folder
    """
        # Convert lists into np.arrays
        frames        = np.arange(0, len(information), 1)
        slb_distances = np.array([entry['slb_dist' ] for entry in information])
        adsorptions   = np.array([entry['site_type' ] for entry in information])
        coordinations = np.array([entry['tot_coord'] for entry in information])
        interactions  = np.array([entry['ads_type' ] for entry in information])


        fig = plt.figure(figsize=(25.08, 17.08), facecolor='w',constrained_layout=True)
        gs  = fig.add_gridspec(2,2,wspace=0.04,hspace=0.05)                                     # divide fig in 4

        # -------- Top-left plot: evolution of distances over time --------
        ax1 = fig.add_subplot(gs[0,0])                                  # upper-left
        point_colors = [site_colors.get(site_color) for site_color in adsorptions]
        ax1.scatter(frames, slb_distances, c = point_colors, alpha = 0.5,
                            linewidth = 0.5)

        ax1.set_xlim(frames[0], frames[-1])
        ax1.set_ylim(0, np.linalg.norm(traj[-1].cell[2])) # assumes cubic cell
        ax1.set_xlabel('MD steps', fontdict=font)
        ax1.set_ylabel('Surface distance (Å)', fontdict=font)

        legend_labels = list(site_colors.keys())
        legend_handles = [plt.Line2D([0], [0], marker='o', color='w', alpha = 0.5,
                                          label=label, markerfacecolor=site_colors[label],
                                          markersize=10) for label in legend_labels]
        ax1.legend(handles=legend_handles, labels=legend_labels, loc='upper left',
                           fontsize=25)

        ax1.set_xlim(frames[0], frames[-1])
        ax1.set_xticks(range(0, frames[-1], int(len(frames)/5)))
        ax1.set_xticklabels(range(0, frames[-1], int(len(frames)/5)))
        for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(30)
        for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(30)

        # -------- Top-right plot: hystogram of distances --------
        ax2 = fig.add_subplot(gs[0,1])
        # Count occurrences of each adsorbant type
        unique_adsorptions, counts = np.unique(adsorptions, return_counts=True)
        # Plot each bar with corresponding color
        for ad, count in zip(unique_adsorptions, counts):
                ax2.bar(ad, count, color=site_colors.get(ad, 'black'), alpha=0.5)

        ax2.set_xlabel('Bonding regime', fontdict=font)
        ax2.set_ylabel('Frequency', fontdict=font)

        for tick in ax2.xaxis.get_major_ticks():
                tick.label1.set_fontsize(30)
        for tick in ax2.yaxis.get_major_ticks():
                tick.label1.set_fontsize(30)

        # -------- Bottom-left plot: evolution of coordinations over time --------
        ax3 = fig.add_subplot(gs[1,0])
        point_colors = [itr_colors.get(ads_color) for ads_color in interactions]
        ax3.scatter(frames, coordinations, c = point_colors, alpha = 0.5,
                            linewidth = 0.5)

        ax3.set_xlabel('MD steps', fontdict=font)
        ax3.set_ylabel('Coordination', fontdict=font)

        legend_labels = list(itr_colors.keys())
        legend_handles = [plt.Line2D([0], [0], marker='o', color='w', alpha = 0.5,
                                          label=label, markerfacecolor=itr_colors[label],
                                          markersize=10) for label in legend_labels]
        ax3.legend(handles=legend_handles, labels=legend_labels, loc='upper left',
                           fontsize=25)

        ax3.set_xlim(frames[0], frames[-1])
        ax1.set_xticks(range(0, frames[-1], int(len(frames)/5)))
        ax1.set_xticklabels(range(0, frames[-1], int(len(frames)/5)))
        ax3.set_yticks([0, 1, 2, 3, 4, 5, 6])
        ax3.set_yticklabels([0, 1, 2, 3, 4, 5, 6])
        for tick in ax3.xaxis.get_major_ticks():
                tick.label1.set_fontsize(30)
        for tick in ax3.yaxis.get_major_ticks():
                tick.label1.set_fontsize(30)

        # -------- Bottom-right plot: hystogram of coordinates --------
        ax4 = fig.add_subplot(gs[1,1])
        unique_interactions, counts = np.unique(interactions, return_counts=True)
        #Plot each bar with corresponding color
        for it, count in zip(unique_interactions, counts):
            ax4.bar(it, count, color=itr_colors.get(it, 'black'), alpha=0.5)

        ax4.set_xlabel('Interaction type', fontdict=font)
        ax4.set_ylabel('Frequency', fontdict=font)

        for tick in ax4.xaxis.get_major_ticks():
                tick.label1.set_fontsize(30)
        for tick in ax4.yaxis.get_major_ticks():
                tick.label1.set_fontsize(30)

        #title
        #plt.suptitle(f'{ads}: atom #{Ai}', fontsize = 30)
        #plt.subplots_adjust(top=0.90)
        plt.savefig(f'Figs/{ads}-{Ai}.png')
        print('Figure saved in the folder Figs!\n')
        plt.close()


# ------------ MAIN --------------- #
# Get user input
traj, frames, ads, ads_indeces, slb, slb_indeces = parse_input()
# all adsorbant will be considered in the connectivity
all_ads_atoms = [i for i, atom in enumerate(frames[-1])
                      if atom.symbol == ads]

# Create two folders where data and figs are gonna be stored
os.system('mkdir Data')
os.system('mkdir Figs')

# Huge list containing all infos of each adsorbant per frame in ROW
all_ads_infos = []


# First loop: run through each frame of the trajectory and
# compute connectivity  matrix
for frame in frames:
        all_ads_infos_per_frame = []
        positions  = frame.get_positions()
        all_coords = get_connectivity_matrix(frame)
        # Second loop: run through the i-th adsorbant index, Ai, individually
        for Ai in ads_indeces:
                # Find all distances of atom Ai w.r.t all others, then look for the
                # minimum value filtering among the slab indeces
                all_distances  = get_distances(positions[Ai], positions,
                                                                       cell=frame.cell, pbc=frame.pbc)[1][0]
                slb_dist = np.min(all_distances[slb_indeces])

                # Find through the connectivity matrix the amount of indeces that are
                # coordinated as molecule and as adsorbant-surface interactions
                coords_indeces = np.nonzero(all_coords[Ai])[0]
                in_ads = set(coords_indeces).intersection(set(all_ads_atoms))
                in_slb = set(coords_indeces).intersection(set(slb_indeces))
                tot_coord = np.sum(all_coords[Ai])
                try:
                        site_type  = ads_sites[str(len(in_slb))]
                except:
                        site_type = 'in bulk'

                if len(in_ads) == 0:
                        mol_type = f'{ads}'
                        if len(in_slb) == 0:
                                ads_type =  'gas H'
                        else:
                                ads_type = 'H*'
                else:
                        mol_type = f'{ads}{len(in_ads)+1}'
                        if len(in_slb) == 0:
                                ads_type = 'gas H\u2099'
                        else:
                                ads_type = 'H\u2099*'



                # Create a dictionary with information about Ai in frame
                Ai_infos = {'slb_dist'  : slb_dist,
                                'tot_coord' : tot_coord,
                    'site_type' : site_type,
                                        'mol_type'  : mol_type, # not used in the plotting
                                        'ads_type'  : ads_type,
                                        }
                all_ads_infos_per_frame.append(Ai_infos)
        all_ads_infos.append(all_ads_infos_per_frame)

# Now each row represents each ads evolution over time
all_ads_infos = np.array(all_ads_infos).T


for i in range(0, len(all_ads_infos)):
        Ai            = ads_indeces[i]
        Ai_ads_infos  = all_ads_infos[i]
        write_json(f'Data/{ads}-{Ai}_infos.json', Ai_ads_infos)
        print(f'Arrays saved in Data/{ads}-{Ai}_infos.dat')

        # Plot time-evolution of ads-(Ai)
        plot_stats(Ai_ads_infos)