# Regression Tests: Python vs Julia TauRunner

This directory contains regression tests comparing the Julia TauRunner.jl implementation against the original Python TauRunner.

## Scripts

| Purpose | Script |
|---------|--------|
| Generate Python reference data | `generate_python_reference.py` |
| Compare implementations | `compare_implementations.jl` |
| Reference data file | `reference_data.json` |

## Running the Tests

### Generate Python Reference Data

Requires Python with `numpy`, `scipy`, and `proposal` installed:

```bash
python generate_python_reference.py --output reference_data.json
```

On ARM Macs with x86_64 PROPOSAL bindings, use Rosetta:

```bash
arch -x86_64 /usr/bin/python3 generate_python_reference.py --output reference_data.json
```

Options:
- `--skip-mc`: Skip Monte Carlo tests (geometry only)
- `--seed N`: Set random seed for MC (default: 42)

### Run Julia Comparison

```bash
julia --project=../.. compare_implementations.jl reference_data.json
```

## Results Summary

| Category | Tests | Max Relative Difference | Status |
|----------|-------|------------------------|--------|
| Earth Radius | 1 | 0% (exact) | PASS |
| Density Profile | 11 | 0% (exact) | PASS |
| Chord Length | 10 | <10^-12 % | PASS |
| Column Depth | 10 | <0.13% | PASS |
| x_to_r Transform | 50 | 0% (exact) | PASS |

## Detailed Results

### Earth Radius

| Parameter | Python | Julia | Rel. Diff |
|-----------|--------|-------|-----------|
| Radius (km) | 6368.0 | 6368.0 | **0%** |

### Density Profile (g/cm^3)

| r | Python | Julia | Rel. Diff |
|---|--------|-------|-----------|
| 0.0 | 13.0885 | 13.0885 | 0% |
| 0.1 | 13.0001 | 13.0001 | 0% |
| 0.2 | 12.1388 | 12.1388 | 0% |
| 0.3 | 11.7253 | 11.7253 | 0% |
| 0.4 | 11.1394 | 11.1394 | 0% |
| 0.5 | 10.3479 | 10.3479 | 0% |
| 0.6 | 5.3956 | 5.3956 | 0% |
| 0.7 | 5.0754 | 5.0754 | 0% |
| 0.8 | 4.7364 | 4.7364 | 0% |
| 0.9 | 3.9845 | 3.9845 | 0% |
| 1.0 | 2.6000 | 2.6000 | 0% |

### Chord Total Length (units of radius)

| theta (deg) | Python | Julia | Rel. Diff |
|-------------|--------|-------|-----------|
| 0 | 2.000000 | 2.000000 | 0% |
| 10 | 1.969616 | 1.969616 | 0% |
| 20 | 1.879385 | 1.879385 | 0% |
| 30 | 1.732051 | 1.732051 | 0% |
| 45 | 1.414214 | 1.414214 | 0% |
| 60 | 1.000000 | 1.000000 | 0% |
| 70 | 0.684040 | 0.684040 | 0% |
| 80 | 0.347296 | 0.347296 | 0% |
| 85 | 0.174311 | 0.174311 | <10^-14 % |
| 89 | 0.034905 | 0.034905 | <10^-12 % |

### Column Depth (g/cm^2)

| theta (deg) | Python | Julia | Rel. Diff |
|-------------|--------|-------|-----------|
| 0 | 1.0945e+10 | 1.0945e+10 | 0.000% |
| 10 | 1.0410e+10 | 1.0409e+10 | 0.010% |
| 20 | 9.0342e+09 | 9.0342e+09 | 0.000% |
| 30 | 6.7907e+09 | 6.7907e+09 | 0.000% |
| 45 | 4.0992e+09 | 4.0992e+09 | 0.000% |
| 60 | 2.5526e+09 | 2.5526e+09 | 0.000% |
| 70 | 1.4971e+09 | 1.4970e+09 | 0.009% |
| 80 | 7.2983e+08 | 7.2983e+08 | 0.000% |
| 85 | 3.3089e+08 | 3.3049e+08 | 0.122% |
| 89 | 5.7792e+07 | 5.7791e+07 | 0.001% |

### x_to_r Coordinate Transform

All 50 test points (10 angles x 5 x-values) match exactly (0% relative difference).

## Monte Carlo Reference Data

Reference data from Python TauRunner (seed=42):

| E (eV) | theta (deg) | N | nu_tau survived | mean E_out (eV) | tau exited | mean n_CC | mean n_NC |
|--------|-------------|---|-----------------|-----------------|------------|-----------|-----------|
| 1e+15 | 85 | 100 | 100 | 9.08e+14 | 0 | 0.11 | 0.06 |
| 1e+15 | 60 | 100 | 100 | 4.84e+14 | 0 | 0.63 | 0.23 |
| 1e+15 | 30 | 100 | 100 | 1.31e+14 | 0 | 1.15 | 0.42 |
| 1e+17 | 85 | 100 | 100 | 4.61e+16 | 0 | 0.68 | 0.24 |
| 1e+17 | 60 | 100 | 100 | 1.75e+15 | 0 | 1.62 | 0.56 |
| 1e+18 | 89 | 100 | 98 | 6.78e+17 | 2 | 0.31 | 0.16 |
| 1e+18 | 85 | 100 | 100 | 1.75e+17 | 0 | 0.91 | 0.31 |

**Note:** Full MC comparison requires cross-sections to be loaded in Julia. The reference data is captured for future validation when cross-section support is complete.

## Precision Notes

| Quantity | Agreement | Notes |
|----------|-----------|-------|
| Geometry (radius, density, coordinates) | Exact | Floating-point identical |
| Chord lengths | <10^-12 % | Machine precision |
| Column depths | <0.13% | Numerical integration tolerance |

The worst-case column depth difference (0.12% at theta=85 deg) is due to different numerical integration methods:
- **Python:** `scipy.integrate.quad` with cached 2D spline interpolation
- **Julia:** `QuadGK.quadgk` with direct integration

Both methods are well within acceptable tolerance for physics applications.
