import numpy as np
from scipy.special import sph_harm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import ase.io
import os

class Featurizer:
    def __init__(self, r_max=8.0, n_radial=8, max_ell=3, cutoff_power=5):
        """
        Feature calculator with rotational invariance

        Parameters:
            r_max: float (default=8.0)
                Maximum neighbor distance
            n_radial: int (default=8)
                Number of radial basis functions
            max_ell: int (default=3)
                Maximum spherical harmonic degree
            cutoff_power: int (default=5)
                Power for polynomial cutoff function
        """
        self.r_max = r_max
        self.n_radial = n_radial
        self.max_ell = max_ell
        self.cutoff_power = cutoff_power

        # Bessel wavenumbers for l=0
        self.bnu0 = np.arange(1, n_radial + 1) * np.pi / r_max

        # number of features
        # n_radial from l=0 and nradial*(nradial+1)/2 from each (max_ell of them) l>0
        self.feature_dim = self.n_radial + self.max_ell * (self.n_radial * (self.n_radial + 1)) // 2

        # Indices for upper triangles after contraction over m
        self.upper_triangle = np.triu_indices(self.n_radial)

    def compute_features(self, atoms):
        """Computes features on ASE atoms object"""
        features = np.zeros((len(atoms), self.feature_dim))
        pos = atoms.get_scaled_positions()
        cell = np.array(atoms.cell.tolist())

        for i,p in enumerate(pos):
            # minumum image convention interatomic distances
            dist = pos - p
            dist -= np.rint(dist)
            rij = np.dot(dist, cell)
            rad = np.linalg.norm(rij, axis = 1)
            neighb = rad <= self.r_max
            neighb[i] = False
            r = rad[neighb]

            if len(r) == 0: # lonely atom
                continue

            # Radial basis and cutoff
            radial = np.sin(np.multiply.outer(self.bnu0, r)) / r \
                   * np.sqrt(2 / self.r_max) * self._cutoff(r)

            # Get spherical angles
            theta, phi = self._cartesian_to_spherical(rij[neighb])

            # Compute density terms
            rho = []
            for l in range(self.max_ell + 1):
                Ylm = self._real_spherical_harmonics(l, theta, phi)
                rho.append(radial @ Ylm.T)

            # Compute invariants
            features[i] = self._compute_invariants(rho)

        return features

    def _cutoff(self, r):
        """Polynomial cutoff, adapted from MACE's radial.py"""
        x = r / self.r_max
        p = self.cutoff_power
        return (1.0 - ((p + 1.0) * (p + 2.0) / 2.0) * x**p \
               + p * (p + 2.0) * x**(p+1) \
               - (p * (p + 1.0) / 2) * x**(p+2)) * (r < self.r_max)

    def _cartesian_to_spherical(self, vectors):
        """Convert Cartesian to spherical coordinates"""
        x, y, z = vectors.T
        r = np.linalg.norm(vectors, axis=1)
        theta = np.arccos(np.divide(z, r, out=np.zeros_like(z), where=r>1e-12))
        phi = np.arctan2(y, x) % (2*np.pi)
        return theta, phi

    def _real_spherical_harmonics(self, l, theta, phi):
        """Compute real spherical harmonics"""
        Ylm = []
        for m in range(-l, 0):
            y = sph_harm(-m, l, phi, theta)
            Ylm.append(np.sqrt(2) * (-1)**m * y.imag)
        y = sph_harm(0, l, phi, theta)
        Ylm.append(y.real)
        for m in range(1, l+1):
            y = sph_harm(m, l, phi, theta)
            Ylm.append(np.sqrt(2) * (-1)**m * y.real)
        return np.stack(Ylm, axis=0)

    def _compute_invariants(self, rho):
        """Compute invariant features"""
        # l=0:
        features = [rho[0][:, 0]]

        # contract over m for l>0
        for l in range(1, self.max_ell + 1):
            features.append((rho[l] @ rho[l].T)[self.upper_triangle])

        return np.concatenate(features)


def build_environment_matrix(features, species, n_clusters=100):
    """Construct matrix C by clustering embeddings."""
    # Flatten all embeddings and species
    all_feat = np.vstack(features)
    all_z = np.concatenate(species)
    unique_z = np.unique(all_z)

    # Cluster per atomic species
    scalers, kmeans = {}, {}
    label_offset = 0
    env_types = {}

    for z in unique_z:
        mask = (all_z == z)
        X = all_feat[mask]

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Cluster
        km = KMeans(n_clusters=n_clusters, n_init=5, random_state=123456)
        km.fit(X_scaled)

        scalers[z] = scaler
        kmeans[z] = km
        env_types[z] = list(range(label_offset, label_offset + n_clusters))
        label_offset += n_clusters

    # Assign environment labels per structure
    C = np.zeros((len(features), label_offset), dtype=int)

    for i, (ft, z) in enumerate(zip(features, species)):
        labels_z = np.zeros(len(z), dtype=int)
        for atom_z in unique_z:
            mask = (z == atom_z)
            if not np.any(mask):
                continue
            X_scaled = scalers[atom_z].transform(ft[mask])
            labels_z[mask] = kmeans[atom_z].predict(X_scaled) + env_types[atom_z][0]
        C[i, :] = np.bincount(labels_z, minlength=label_offset)

    return C



if __name__ == "__main__":
    # Cache file for environmen matrix
    matrix_file = 'environment_matrix.npz'

    # Check if the file exists
    if os.path.exists(matrix_file):
        # Load the matrix from file
        print(f"Loading matrix from {matrix_file}")
        data = np.load(matrix_file)
        C = data['C']
    else:
        # Compute the matrix if it doesn't exist
        print(f"Computing matrix and saving to {matrix_file}")
        feat = Featurizer(r_max=8.0, n_radial=8, max_ell=3, cutoff_power=5)

        # Read training & testing configurations
        atomsdata = ase.io.read('struct.xyz', index=':')

        # Compute features
        features = []
        species = []
        for atoms in atomsdata:
            features.append(feat.compute_features(atoms))
            species.append(atoms.get_atomic_numbers())

        # Compute environment matrix
        C = build_environment_matrix(features, species, n_clusters=100)

        # Save the matrix
        np.savez(matrix_file, C=C)

    # Compute SVD and rank of C.
    U, s, Vt = np.linalg.svd(C, full_matrices=False)
    rank = np.linalg.matrix_rank(C)

    print(f"Matrix shape: {C.shape}")
    print(f"Numerical rank: {rank}")
    print(f"Top 10 singular values:\n{s[:10]}")
    print(f"Lowest 10 singular values:\n{s[-10:]}")

    # Plot singular values
    plt.figure(figsize=(8, 4))
    plt.plot(s, 'o-', markersize=4)
    plt.axvline(rank, color='r', linestyle='--', label=f'Rank = {rank}')
    plt.xlabel('Index')
    plt.ylabel('Singular Value')
    plt.yscale('log')
    plt.title('Singular Value Spectrum of Environment Matrix C')
    plt.legend()
    plt.tight_layout()
    plt.savefig('singular_values.png')
    plt.close()

    # Environment matrix plot
    plt.style.use('seaborn-v0_8-notebook')

    fig, ax = plt.subplots(figsize=(12, 8), dpi=120)

    # Colormap
    cmap = plt.cm.plasma

    # Plot matrix
    im = ax.imshow(C.T,
                  cmap=cmap,
                  origin='lower',
                  aspect='auto',
                  interpolation='nearest')

    # Add colorbar for environment counts
    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label('Environment Counts',
                  fontsize=20,
                  rotation=270,
                  labelpad=30)
    cbar.ax.tick_params(labelsize=20)

    ax.set_xlabel('Configuration Index',
                 fontsize=20, 
                 labelpad=10)
    ax.set_ylabel('Environment Cluster', 
                 fontsize=20, 
                 labelpad=10)
    ax.set_title('Environment Count Matrix\n(Per-configuration cluster counts)',
               fontsize=20,
               pad=20)

    ax.tick_params(axis='both', which='both', labelsize=20)
    ax.grid(False)  # Turn off grid for matrix plot
    ax.minorticks_on()
    plt.tight_layout()

    plt.savefig('environment_matrix.pdf', bbox_inches='tight')
