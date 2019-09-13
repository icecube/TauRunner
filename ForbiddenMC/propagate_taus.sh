#!/bin/bash

Array=( "$@" )
last_idx=$(( ${#Array[@]} - 1 ))
rho=${Array[$last_idx]}
unset Array[$last_idx]

DIR='/data/user/apizzuto/ANITA/monte_carlo/code_cleanup/TauDragon/ForbiddenMC/'

awk -v var="${Array[*]}" "BEGIN {n=split(var , array); for (i=1; i<=n; i++) print array[i], 100000000}" | ${DIR}../MMC/MMC/ammc -run -frejus -tau -medi="Frejus rock" -radius=1e6 -vcut=1.e-3 -rho=$rho -scat -lpm -bs=1 -ph=3 -bb=1 -sh=1 -ebig=1e16 -seed=1223 -tdir=${DIR}../MMC/tables/

