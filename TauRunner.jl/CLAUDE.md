# TauRunner.jl Development Notes

## Project Overview
TauRunner.jl is a Julia port of the Python TauRunner package for Monte Carlo simulation of ultra-high-energy neutrino propagation through matter (Earth, Sun, arbitrary spherical bodies).

## Dependencies

### PROPOSAL.jl (Local Development Package)
The PROPOSAL library bindings for Julia are located at:
```
~/research/proposal_julia/PROPOSAL.jl
```

This package is registered as a local development dependency. If you need to re-add it:
```julia
using Pkg
Pkg.develop(path="/Users/jlazar/research/proposal_julia/PROPOSAL.jl")
```

The upstream PROPOSAL C++ library is located at:
```
~/research/proposal_julia/PROPOSAL_upstream/
```

### PROPOSAL Interpolation Tables
PROPOSAL generates interpolation tables on first use. By default, TauRunner.jl configures tables to be stored in:
```
data/proposal_tables/
```

You can override this by setting the `PROPOSAL_TABLES_PATH` environment variable:
```bash
export PROPOSAL_TABLES_PATH=/path/to/your/tables
julia --project=. scripts/compare_cross_sections.jl
```

**Note:** The JSON config `tables_path` setting may not be respected by all PROPOSAL versions. If tables are still written to `/tmp`, the PROPOSAL.jl bindings may need to expose `InterpolationSettings.tables_path` directly (like Python PROPOSAL does).

### Known Differences from Python TauRunner
The differential cross-section splines (`dσ/dz`) are interpolated using **bilinear** interpolation in Julia (`Interpolations.jl linear_interpolation`) versus **bicubic** splines in Python (`scipy.interpolate.RectBivariateSpline`, which defaults to `kx=3, ky=3`). Both use the same underlying 200x200 grid of knot points (stored as npz/pickle). The bilinear interpolation produces a slightly higher mean inelasticity `z` (i.e., the neutrino retains more energy per interaction), with the effect growing at higher energies:

| log10(E/eV) | Δmean(z) bilinear−bicubic |
|-------------|---------------------------|
| 15          | 0.0003                    |
| 16          | 0.003                     |
| 17          | 0.012                     |
| 18          | 0.029                     |

This compounds over multiple interactions, producing a ~0.01−0.02 dex systematic shift toward higher outgoing energies in Julia relative to Python. The total cross-sections, column depth calculations, interaction sampling logic, and PROPOSAL charged lepton propagation all agree between the two implementations.

### Known Issues
- The PROPOSAL spherical body propagator requires an outer "Air" sector to handle particles starting at the surface. This has been fixed in `ConfigGeneration.jl`.
- PROPOSAL may show numerical warnings ("Precision not reached", "Maximum iteration exceeded in Bisection") which are generally harmless.

## Running Tests
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Comparison Scripts
To compare Julia output with Python TauRunner:
```bash
julia --project=. scripts/compare_cross_sections.jl
```

This script compares:
- Units and conversion factors
- Physics constants
- Bodies/Geometry (Earth radius, density profile)
- Tracks (Chord properties, column depth)
- Particle types
- Cross-sections (total and differential, CSMS and DIPOLE models)
- Monte Carlo propagation (interaction depths, single particle propagation, statistics)

## Directory Structure
- `src/` - Main source code
- `data/cross_sections/` - Cross-section data files (JLD2 format)
- `data/proposal_tables/` - PROPOSAL interpolation tables (generated on first use)
- `scripts/` - Utility and comparison scripts
- `test/` - Test suite including regression tests
- `examples/` - Usage examples
