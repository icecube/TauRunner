# TauRunner.jl

A Julia implementation of TauRunner for propagating ultra-high-energy neutrinos through matter.

## Overview

TauRunner.jl is a Julia port of the [TauRunner](https://github.com/icecube/TauRunner) Python package. It propagates neutrinos (with a focus on tau neutrinos) through spherical bodies like Earth using Monte Carlo methods, accounting for:

- Neutrino-nucleon interactions (CC and NC)
- Tau lepton energy losses and decay
- Regeneration from tau decay
- PREM Earth density model

## Installation

Requires Julia 1.11 or later.

PROPOSAL.jl is not yet registered in the General registry, so it must be added explicitly before TauRunner.jl:

```julia
using Pkg
Pkg.add(url="https://github.com/jlazar17/PROPOSAL.jl")
Pkg.add(url="https://github.com/icecube/TauRunner", subdir="TauRunner.jl")
```

For development:
```julia
using Pkg
Pkg.develop(path="path/to/TauRunner/TauRunner.jl")
```

## Quick Start

```julia
using TauRunner

# Create Earth with PREM density model
earth = construct_earth()

# Create a chord through the Earth (theta=0 goes through the core)
chord = Chord(0.0)  # nadir angle in radians

# Compute column depth
depth = total_column_depth(chord, earth)
println("Column depth: ", depth / (TauRunner.units.gr / TauRunner.units.cm^2), " g/cm²")
```

## Features

- **Bodies**: PREM Earth model, custom spherical bodies, slab geometries
- **Tracks**: Chord trajectories, radial trajectories, slab tracks
- **Cross-sections**: DIPOLE and CSMS neutrino-nucleon cross-section models
- **Propagation**: Full Monte Carlo with interaction sampling and tau regeneration
- **Caching**: LRU-cached column depth calculations for performance

## Testing

```julia
using Pkg
Pkg.test("TauRunner")
```

## Known Differences from Python TauRunner

The differential cross-section splines (dσ/dz) are interpolated using bilinear interpolation in Julia (`Interpolations.jl`) versus bicubic splines in Python (`scipy.interpolate.RectBivariateSpline`). Both use the same underlying 200x200 grid of knot points. The bilinear interpolation produces a slightly higher mean inelasticity z (the neutrino retains more energy per interaction), with the effect growing at higher energies:

| log10(E/eV) | Δmean(z) bilinear−bicubic |
|-------------|---------------------------|
| 15          | 0.0003                    |
| 16          | 0.003                     |
| 17          | 0.012                     |
| 18          | 0.029                     |

This compounds over multiple interactions, producing a ~0.01−0.02 dex systematic shift toward higher outgoing energies in Julia relative to Python. The total cross-sections, column depth calculations, interaction sampling logic, and PROPOSAL charged lepton propagation all agree between the two implementations.

## Authors

Ibrahim Safa, Carlos A. Arguelles, Jeff Lazar, Alex Pizzuto

## Citation

If you use TauRunner in your research, please cite:

> Ibrahim Safa, Jeffrey Lazar, Alex Pizzuto, Oswaldo Vasquez, Carlos A. Argüelles, Justin Vandenbroucke
> "TauRunner: A Public Python Program to Propagate Neutral and Charged Leptons"
> Comput. Phys. Commun. 278 (2022) 108422. [arXiv:2110.14662](https://arxiv.org/abs/12110.14662)
