# TauRunner.jl

A Julia implementation of TauRunner for propagating ultra-high-energy neutrinos through matter.

## Overview

TauRunner.jl is a Julia port of the [TauRunner](https://github.com/icecube/TauRunner) Python package. It propagates neutrinos (with a focus on tau neutrinos) through spherical bodies like Earth using Monte Carlo methods, accounting for:

- Neutrino-nucleon interactions (CC and NC)
- Tau lepton energy losses and decay
- Regeneration from tau decay
- PREM Earth density model

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/icecube/TauRunner", subdir="TauRunner.jl")
```

Or for development:
```julia
using Pkg
Pkg.develop(path="path/to/TauRunner/TauRunner.jl")
```

## PROPOSAL.jl Dependency

TauRunner.jl requires [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL) for charged lepton propagation. The Julia bindings (PROPOSAL.jl) are currently **not yet available** through Yggdrasil, but the submission is pending. Until then, PROPOSAL.jl must be installed as a local development dependency.

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

> Ibrahim Safa, Alex Pizzuto, Carlos Arguelles, Francis Halzen, Raamis Hussain, Ali Kheirandish, Justin Vandenbroucke.
> "Observing EeV neutrinos through Earth: GZK and the anomalous ANITA events"
> JCAP 01 (2020) 012. [arXiv:1909.10487](https://arxiv.org/abs/1909.10487)
