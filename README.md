# StagCorr

A Python package for computing staggered fermion correlation functions in lattice QCD simulations.

## Overview

StagCorr provides tools for calculating n-point correlation functions using the staggered fermion formalism in lattice Quantum Chromodynamics (QCD). It supports both free field theory calculations and full QCD simulations with gauge field interactions, with results validated against MILC reference data.

## Features

- **N-point correlation functions**: 2-point, 3-point, and 4-point correlators
- **Staggered fermion formalism**: Complete implementation of staggered fermion phases and operators  
- **Gauge field support**: HISQ smearing and gauge field handling
- **Flexible lattice geometries**: Support for various lattice sizes and boundary conditions
- **Free field theory**: Analytical calculations for comparison and testing
- **MILC validation**: Results validated against MILC lattice QCD code

## Installation

### From PyPI (when published)
```bash
pip install stagcorr
```

### From source
```bash
git clone https://github.com/shaun252/stagcorr.git
cd stagcorr
pip install -e .
```

## Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20.0
- SymPy ≥ 1.8.0

## Quick Start

```python
import stagcorr.correlators.main as corr

# Set lattice parameters
vol = (6, 6, 6, 6)  # (T,X,Y,Z) dimensions

# Generate a pion two-point correlation function
Ct = corr.generate_npt(
    spinTasteMassNaikMomSymShift1=["G5", "G5", 1, 0, [0,0,0], 0],
    spinTasteMassNaikMomSymShift2=["G5", "G5", 1, 0, [0,0,0], 0],
    volume=vol
)

print(Ct[:, 0])  # Print correlation function values
```

## Examples

### Two-Point Functions

```python
# Vector meson correlator: GX⊗G1 × GX⊗G1
# Mass = 1.0 (heavy quark), Naik epsilon = 0 (no improvement)
Ct = corr.generate_npt(
    spinTasteMassNaikMomSymShift1=["GX", "G1", 1.0, 0, [0,0,0], 0],
    spinTasteMassNaikMomSymShift2=["GX", "G1", 1.0, 0, [0,0,0], 0],
    volume=vol
)
# Returns: nt×nt tensor indexed as Ct[t2, t1]
```

### Three-Point Functions

```python
# Two-pion to vector transition: π(+p) + π(-p) → ρ(0)
# Mass = 0.1 (light quark), Naik epsilon = 0 (no improvement)
Ct = corr.generate_npt(
    spinTasteMassNaikMomSymShift1=["G5", "G5", 0.1, 0, [0,0,1], 0],   # π⁺
    spinTasteMassNaikMomSymShift2=["G5", "G5", 0.1, 0, [0,0,-1], 0],  # π⁻ 
    spinTasteMassNaikMomSymShift3=["GZ", "G1", 0.1, 0, [0,0,0], 0],   # ρ
    volume=vol
)
# Returns: nt×nt×nt tensor indexed as Ct[t3, t2, t1]
```

### Four-Point Functions

```python
# Two-pion scattering: π(+p) + π(-p) → π(+p) + π(-p)
# Mass = 0.1 (light quark), Naik epsilon = 0 (no improvement)
Ct = corr.generate_npt(
    spinTasteMassNaikMomSymShift1=["G5", "G5", 0.1, 0, [0,0,1], 0],   # incoming π⁺
    spinTasteMassNaikMomSymShift2=["G5", "G5", 0.1, 0, [0,0,-1], 0],  # incoming π⁻
    spinTasteMassNaikMomSymShift3=["G5", "G5", 0.1, 0, [0,0,1], 0],   # outgoing π⁺
    spinTasteMassNaikMomSymShift4=["G5", "G5", 0.1, 0, [0,0,-1], 0],  # outgoing π⁻
    volume=vol
)
# Returns: nt×nt×nt×nt tensor indexed as Ct[t4, t3, t2, t1]
```

## Parameter Format

Each operator is specified as `[spin, taste, mass, naik_epsilon, momentum, symmetric_shift]`:

- **Spin**: Dirac structure (`"G5"`, `"GX"`, `"GY"`, `"GZ"`, `"GT"`, `"G5X"`, etc.)
- **Taste**: Taste structure (`"G5"`, `"G1"`, `"GX"`, `"GY"`, `"GZ"`, etc.)
- **Mass**: Quark mass in lattice units (0.1 = light, 1.0 = heavy)
- **Naik epsilon**: Naik improvement parameter (0 = no improvement, typically 0)
- **Momentum**: `[px, py, pz]` in lattice units (2π/L)
- **Symmetric shift**: Shift parameter (typically 0)

## Output Format

The `generate_npt()` function returns correlation function tensors:

- **2-point**: `nt×nt` tensor indexed as `Ct[t2, t1]`
- **3-point**: `nt×nt×nt` tensor indexed as `Ct[t3, t2, t1]`  
- **4-point**: `nt×nt×nt×nt` tensor indexed as `Ct[t4, t3, t2, t1]`

where `nt` is the temporal lattice size and times are ordered chronologically (t1 < t2 < t3 < t4).

## Operator Notation

Following staggered fermion conventions:
- `G5` = γ₅, `GX` = γ₁, `GY` = γ₂, `GZ` = γ₃, `GT` = γ₀
- `G5X` = γ₅γ₁, `GXT` = γ₁γ₀, etc.
- `G1` = Identity (taste singlet)

## Package Structure

- `stagcorr.correlators`: Core correlation function calculations
  - `main.py`: Primary interface via `generate_npt()` function
  - `twoPoints.py`, `twoPointsFree.py`: 2-point correlators
  - `threePoints.py`, `threePointsFree.py`: 3-point correlators  
  - `fourPoints.py`, `fourPointsFree.py`: 4-point correlators
- `stagcorr.gauge`: Gauge field operations and HISQ smearing
- `stagcorr.lattice`: Lattice navigation utilities
- `stagcorr.propagatorUtils`: Staggered Dirac operator and propagator tools
- `stagcorr.stagFuncs.py`: Staggered fermion phase functions

## Physics Background

This package implements the staggered fermion discretization of QCD on a 4-dimensional Euclidean lattice. Staggered fermions reduce the number of fermion doublers while maintaining chiral symmetry properties, making them particularly useful for lattice QCD calculations.

The package uses:
- 2×stag_op (=1/2 × prop) convention consistent with MILC
- HISQ staggered operators with numpy.linalg.inv for propagator inversion
- Wick contractions for 2-, 3- and 4-point functions

## Testing

Run the test suite using pytest:

```bash
pytest tests/
```

The tests validate calculations against MILC reference data with agreement typically better than 1e-8 relative precision.

## Validation

Results are validated against the MILC lattice QCD code with agreement typically better than 1e-8 relative precision.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

MIT License

## Citation

If you use this package in your research, please cite:
```
Lahert, S. (2024). StagCorr: Staggered fermion correlation function solver for lattice QCD.
```

## Contact

Shaun Lahert - shaun.lahert@gmail.com