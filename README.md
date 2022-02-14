<img src="https://icecube.wisc.edu/~isafa/TauRunner.jpg" alt="logo"
	title="taurunner logo" width="350" height="250" />
 

Authors: Ibrahim Safa, Carlos A. Argüelles, Jeff Lazar, Alex Pizzuto
=======

## Introduction

`TauRunner` is a tool that propagates ultra-high-energy neutrinos, with a focus on tau neutrinos. Although it was developed for extremely high energy (EeV+) applications, it is able to propagate neutrinos from 1 to 10^16 GeV. Note that oscillations are not taken into account at the lowest energies, but they become negligible above 1 TeV.   

## Installation

`taurunner` is available via PyPI. If you are interested in using `taurunner` but do not plan on contributing to the code base, you can install it with
```console
pip install taurunner
```

If you would like to be able to modify the code, then we recommend downloading the code via git and then installing with `pip`
```console
git clone git@github.com:icecube/TauRunner.git
pip install taurunner
```

If you find any problems with the code or have any feature requests, please open an issue on [Github](https://github.com/icecube/TauRunner/issues).

## Examples
Below are three prototypical command-line use cases for this project. For more examples, see `notebooks/TauRunner_examples.ipynb`.

### Injecting a beam of monoenergetic tau-neutrinos
We can simulate 500 tau neutrinos with an energy of 1 EeV that travel directly through the core. At the point of emergence from the Earth, we receive the resultant particle ID and energy.
```console
python taurunner/main.py -n 500 -t 0.0 -e 1e18
```
`-n` specifies the number of events, `-t` is the nadir angle in degrees (0-90), `-s` specifies a unique seed, and `-e` is the particle energy in eV. 

### Injecting a beam of multi-energy tau-neutrinos following a power-law distribution

We can simulate 500 tau neutrinos with energies following a power-law distribution from 100 Tev to 10 PeV that travel directly through the core. At the point of emergence from the Earth, we receive the resultant particle ID and energy.
```console
python taurunner/main.py -e=-2. -t=30. -n=100 --e_min=1e14 --e_max=1e16
```
`-n` specifies the number of events, `-t` is the nadir angle in degrees, a negative `-e` indicates a power law spectrum with a given spectral index, and `e_min` and `e_max` specify the range over which to sample these energies.

### Simulating a model of an isotropic cosmogenic flux
The user also has the option to simulate an isotropic (over half of the sky) flux, where the energy spectrum is specified by spline of a cdf:
```console
python main.py -e=resources/gzk_cdf_phi_spline.npy -n=10 -t="range" --th_max=90. --th_min=0.
```
if the `-e` argument points to a file, it loads this in and attempts to sample from it assuming the file contains a spline of the CDF of this flux (details on how to construct this in the `notebooks` directory). This is usefull for propagating fluxes with more complicated spectral shapes. 

### Output
Output can either be saved or printed. If saved, data is stored in a `numpy.ndarray`, otherwise it is printed in a format such as the one below (this is for a gzk flux simulated with 5 events):
```console
╒═════════════╤═════════════╤═════════╤═══════╤═══════╤════════════════╕
│        Eini │        Eout │   Theta │   nCC │   nNC │   PDG_Encoding │
╞═════════════╪═════════════╪═════════╪═══════╪═══════╪════════════════╡
│ 8.72753e+17 │ 1.9929e+15  │ 64.6839 │     4 │     1 │             16 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 1.14611e+16 │ 5.29122e+13 │ 67.4882 │     2 │     0 │             16 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 5.53607e+15 │ 2.19641e+13 │ 33.8733 │     3 │     1 │             16 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 4.92978e+17 │ 1.75953e+15 │ 75.4485 │     1 │     0 │             16 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 1.2461e+15  │ 6.1269e+13  │ 67.1521 │     1 │     0 │             16 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 7.76704e+14 │ 7.76704e+14 │ 64.6839 │     0 │     0 │            -14 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 1.23876e+15 │ 0           │ 64.6839 │     1 │     0 │            -11 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 6.85498e+14 │ 0           │ 67.4882 │     1 │     0 │            -13 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 1.43932e+14 │ 1.43932e+14 │ 33.8733 │     0 │     0 │            -14 │
├─────────────┼─────────────┼─────────┼───────┼───────┼────────────────┤
│ 1.29186e+14 │ 1.29186e+14 │ 33.8733 │     0 │     0 │            -14 │
╘═════════════╧═════════════╧═════════╧═══════╧═══════╧════════════════╛
```
The first column provides you with the initial energy, the second with the outgoing energy, the third with the sampled nadir angle. The other columns contain the number of charged current (CC) and neutral current (NC) interactions the particle underwent, and the outgoing particle ID is in the last column.

## Other options
For a full list of the available options available when running `taurunner` from the command line, run
```console
python taurunner/main.py --h
```

## Citation

To cite this work, and for more information, please refer to

Observing EeV neutrinos through Earth: GZK and the anomalous ANITA events

Ibrahim Safa, Alex Pizzuto, Carlos Argüelles, Francis Halzen, Raamis Hussain, Ali Kheirandish, Justin Vandenbroucke

Journal article: https://doi.org/10.1088/1475-7516/2020/01/012
preprint: arXiv:1909.10487

