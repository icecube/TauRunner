<img src="https://icecube.wisc.edu/~isafa/TauRunner.jpg" alt="logo"
	title="taurunner logo" width="350" height="250" />
 

Authors: Ibrahim Safa, Carlos A. Argüelles, Alex Pizzuto
=======

## Introduction

`TauRunner` is a tool that propagates tau-neutrinos and taus through the Earth. Although it was developed for extremely high energy (EeV+) applications, it is able to propagate neutrinos from 1 to 10^16 GeV. Note that oscillations are not taken into account at the lowest energies, but they become negligible above 1TeV.   

## Installation

`TauRunner` has a `nuSQuIDS` dependcy. To run it, you need to: 

1) Install `nuSQuIDS` with python interface. 
  *For more information on `nuSQuIDS` including installation instructions see: https://github.com/arguelles/nuSQuIDS

2) Run `TauRunner/MMC/MMC/ammc -compile -f` to compile the tau propagation code. 
  *For more information on MMC see: https://arxiv.org/abs/hep-ph/0407075

3) You're ready to propagate neutrinos! See the examples section for execution instructions.

## Examples
Below are three prototypical use cases for this project:

### Injecting a beam of monoenergetic tau-neutrinos
We can simulate 500 tau neutrinos with an energy of 1 EeV that travel directly through the core. At the point of emergence from the Earth, we receive the resultant particle ID and energy.
```console
python main.py -n 500 -t 0.0 -s 1 -e 1000000000
```
`-n` specifies the number of events, `-t` is the nadir angle in degrees (0-90), `-s` specifies a unique seed, and `-e` is the particle energy in GeV. 

### Injecting a beam of multi-energy tau-neutrinos following a power-law distribution

We can simulate 500 tau neutrinos with energies following a power-law distribution from 100 Tev to 10 PeV that travel directly through the core. At the point of emergence from the Earth, we receive the resultant particle ID and energy.
```console
python main.py -n 500 -t 0.0 -s 1 -spectrum -2 --range 1e5 1e7
```
`-n` specifies the number of events, `-t` is the nadir angle in degrees, `-s` specifies a unique seed for purposes of reproducibility, `-spectrum` specifies the power-law index, and `--range` gives the range of energies to sample from in GeV. 

### Simulating a model of an isotropic cosmogenic flux
The user also has the option to simulate an isotropic (over half of the sky) flux, where the energy spectrum is specified by spline of a cdf:
```console
python main.py -n 500 -gzk ./gzk_cdf_phi_spline.npy -s 2
```
here, when the `-gzk` option is raised with the path to a cdf spline, the energy spectrum is sampled from the relevant file, and the angle is sampled uniformly over solid angle. An example on how to make one of these cdf splines is shown in `resources/cdf_maker.ipynb`

### Output
Output can either be saved or printed. If saved, data is stored in a `numpy.ndarray`, otherwise it is printed in a format such as the one below (this is for a gzk flux simulated with 5 events):
```console
╒═════════════╤═════════════╤══════════╤═════════════╕
│        Eini │        Eout │    Theta │   CDF_index │
╞═════════════╪═════════════╪══════════╪═════════════╡
│ 2.86713e+16 │ 2.86713e+16 │ 1.26366  │   0.396767  │
├─────────────┼─────────────┼──────────┼─────────────┤
│ 2.01731e+16 │ 2.01731e+16 │ 1.57068  │   0.345561  │
├─────────────┼─────────────┼──────────┼─────────────┤
│ 3.09344e+15 │ 3.85069e+14 │ 1.14063  │   0.0923386 │
├─────────────┼─────────────┼──────────┼─────────────┤
│ 7.88035e+16 │ 4.51484e+13 │ 1.42351  │   0.538817  │
├─────────────┼─────────────┼──────────┼─────────────┤
│ 6.72391e+15 │ 1.44919e+13 │ 0.766526 │   0.18626   │
╘═════════════╧═════════════╧══════════╧═════════════╛
╒════════╤════════╤═════════╤═════════════╕
│ Eini   │ Eout   │ Theta   │ CDF_index   │
╞════════╪════════╪═════════╪═════════════╡
╘════════╧════════╧═════════╧═════════════╛
```
The first column provides you with the initial energy, the second with the outgoing energy, the third with the sampled nadir angle, and finally the cdf index so one can check against the initial energies to make sure their spline is behaving properly. The first table is the tau neutrinos exiting the Earth, and the second table is the taus which exit the Earth, if any.

## Other options
There are a variety of other options not specified in the examples that the user may specify at the command line:
* `spectrum`: instead of using a monoenergetic beam or a splined energy spectrum, simulate a power law (index provided as argument)
* `--range`: If using a power law, this is to specify the range over which to sample energies
* `-spectrum`: instead of using a monoenergetic beam or a splined energy spectrum, simulate a power law (index provided as argument)
* `-buff`: Stop the simulation a finite distance (in kilometers) below the surface of Earth. This is helpful for calculating fluxes incident upon underground detectors.
* `-p`: Path to run script from another directory (if running the script from another directory, provide the path to the script here).
* `-d`: print debug statements during execution
* `-save`: specify the path to where you would like output saved. If no path is provided, output is formatted into a table and printed
* `-water`: int. Add a water layer to the Earth in km. 0 by default.
* `-xs`: string. Choose between two cross section models "dipole" (default) and "CSMS". details of these cross sections are discussed in the linked publication. Note that changing the neutrino cross-section also changes the tau energy loss model consistently (dipole mode for dipole, and ALLM for CSMS).
* `-onlytau`: If only the tau distribution is needed. When this flag is raised, neutrino distributions are not saved.

## Citation

To cite this work, and for more information, please refer to

Observing EeV neutrinos through Earth: GZK and the anomalous ANITA events

Ibrahim Safa, Alex Pizzuto, Carlos Argüelles, Francis Halzen, Raamis Hussain, Ali Kheirandish, Justin Vandenbroucke

Journal article: https://doi.org/10.1088/1475-7516/2020/01/012
preprint: arXiv:1909.10487

