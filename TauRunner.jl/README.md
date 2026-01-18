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

## PROPOSAL.jl Status

**Note:** TauRunner.jl uses [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL) for accurate charged lepton propagation. The Julia bindings (PROPOSAL.jl) are currently **not yet available** through Yggdrasil, but the submission is pending.

Until PROPOSAL.jl is available, TauRunner.jl will use a simplified propagation model for charged leptons. **This simplified model is NOT suitable for physics analysis** - it uses average continuous energy losses instead of proper stochastic sampling.

Once PROPOSAL.jl is available via Yggdrasil, simply add it to your environment and TauRunner.jl will automatically use it:
```julia
using Pkg
Pkg.add("PROPOSAL")
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

## Authors

Ibrahim Safa, Carlos A. Arguelles, Jeff Lazar, Alex Pizzuto

## Citation

If you use TauRunner in your research, please cite:

> Ibrahim Safa, Alex Pizzuto, Carlos Arguelles, Francis Halzen, Raamis Hussain, Ali Kheirandish, Justin Vandenbroucke.
> "Observing EeV neutrinos through Earth: GZK and the anomalous ANITA events"
> JCAP 01 (2020) 012. [arXiv:1909.10487](https://arxiv.org/abs/1909.10487)
