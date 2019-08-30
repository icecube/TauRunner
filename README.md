# TauRunner

Authors: Ibrahim Safa, Carlos A. Arg&uuelles, Alex Pizzuto, Austin Schneider

## Introduction

`TauRunner` is a tool for the propagation of extremely-high-energy tau-neutrinos and taus through the Earth. 

## Installation


## Examples
Below are two prototypical use cases for this project:
### Injecting a beam of monoenergetic tau-neutrinos
We can simulate 500 tau neutrinos with an energy of 1 EeV that travel directly through the core. At the point of emergence from the Earth, we receive the resultant particle ID and energy
```console
python main.py -n 500 -t 0.0 -s 1 -e 1000000000
```
`-n` specifies the number of events, `-t` is the nadir angle in degrees, `-s` specifies a unique seed, and `-e` is for the particle energy in GeV. 

### Simulating a model of an isotropic cosmogenic flux
The user also has the option to simulate an isotropic (over half of the sky) flux, where the energy spectrum is specified by spline of a cdf:
```console
python main.py -n 500 -gzk ./gzk_cdf_phi_spline.npy -s 2
```
here, when the `-gzk` option is raised with the path to a cdf spline, the energy spectrum is sampled from the relevant file, and the angle is sampled uniformly over solid angle. To learn how to make one of these splines, see the example in `notebooks/cdf_spline.ipynb`

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
* `-buff`: Stop the simulation a finite distance (in kilometers) below the surface of Earth. This is helpful for calculating fluxes incident upon underground detectors.
* `-p`: Path to run script from another directory (rarely used, only recommended when working with different versions of the project)
* `-d`: print debug statements at the end of the execution
* `-save`: specify the path to where you would like output saved. If no path is provided, output is formatted into a table and printed

