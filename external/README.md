# External Dependencies

This directory contains two external repositories used in the MLP project, included as Git submodules. These are custom forks of [MACE](https://github.com/ACEsuit/mace) and [ASE](https://gitlab.com/ase/ase).

## Submodules

### 1. [`mace_node/`](https://github.com/sumanbhasker89/mace_node)
- **Base project**: MACE
- **Upstream version**: 0.3.4
- **Modified file**:
  - `mace/calculators/mace.py`

### 2. [`ase_multiproc_calc/`](https://github.com/sumanbhasker89/ase_multiproc_calc)
- **Base project**: ASE
- **Upstream version**: 3.22.1
- **Modified file**:
  - `ase/calculators/mixing.py`
