#!/usr/bin/env python

import pycondor, argparse, sys, os.path
from glob import glob
import numpy as np

error = '/scratch/apizzuto/ANITA/error'
output = '/scratch/apizzuto/ANITA/output'
log = '/scratch/apizzuto/ANITA/log'
submit = '/scratch/apizzuto/ANITA/submit'

job = pycondor.Job('tauRunner_submit_batch_jobs','/data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/main.py',
            error=error,
            output=output,
            log=log,
            submit=submit,
            getenv=True,
            universe='vanilla',
            verbose=2,
            request_memory=10000,
            extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT']
            )
energies = np.logspace(6., 12., 7)
cos_zen = np.linspace(0., 1., 11)
cos_zen = cos_zen[1:]
counter = 0

#for en in energies:
#    for cz in cos_zen:
#        counter += 1
#        theta = np.arccos(cz)
#        job.add_arg('-t {} -s {} -e {} -n 1000000 -p /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/ -save -savedir /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/'.format(theta, counter, en) )


for j in range(100):
    for gzk_model in ['', ' -murase', ' -alosio']:
        counter += 1
        job.add_arg('-gzk{} -s {} -n 1000000 -p /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/ -save -savedir /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/ -buff 10'.format(gzk_model, counter))

#for en in energies:
#    for cz in cos_zen:
#        counter += 1
#        theta = np.arccos(cz)
#        job.add_arg('-t {} -s {} -e {} -n 1000000 -p /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/ -save -savedir /data/user/apizzuto/ANITA/monte_carlo/TauDragon/ForbiddenMC/'.format(theta, counter, en) )


dagman = pycondor.Dagman('batchjobs_tauRunner', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
