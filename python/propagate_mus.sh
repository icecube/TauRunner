#!/bin/bash

Array=( "$@" )
last_idx=$(( ${#Array[@]} - 1 ))
rho=${Array[$last_idx]}
unset Array[$last_idx]

DIR=$PWD/

awk -v var="${Array[*]}" "BEGIN {n=split(var , array); for (i=1; i<=n/2; i++) print array[2*i-1], array[2*i]}" | ${DIR}/../MMC/MMC/ammc -run -frejus -mu -medi="Frejus rock" -radius=1e6 -vcut=1.e-3 -rho=$rho -scat -lpm -bs=1 -ph=3 -bb=1 -sh=1 -ebig=1e21 -tdir=${DIR}../MMC/tables/