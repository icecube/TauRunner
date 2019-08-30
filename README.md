# TauRunner

Authors: Ibrahim Safa, Carlos A. Arg\"uelles, Alex Pizzuto, Austin Schneider

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


## Other features


