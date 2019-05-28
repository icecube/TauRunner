#!/usr/bin/env python

import pycondor, argparse, sys, os.path
from glob import glob
import numpy as np

error = '/scratch/apizzuto/ANITA/error'
output = '/scratch/apizzuto/ANITA/output'
log = '/scratch/apizzuto/ANITA/log'
submit = '/scratch/apizzuto/ANITA/submit'

job = pycondor.Job('cosmogenic_flux_taurunner','/data/user/apizzuto/ANITA/TauDragon/ForbiddenMC/Casino.py',
            error=error,
            output=output,
            log=log,
            submit=submit,
            getenv=True,
            universe='vanilla',
            verbose=2,
            request_memory=3500,
            extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT']
            )

energies = np.logspace(6., 12., 31)
cos_zenith = np.linspace(0., 1., 20)
counter = 0
for en in energies:
    for cz in cos_zenith:
        zenith = np.arccos(cz)
        counter += 1
        job.add_arg('-t {} -s {} -e {} -n 1000000'.format(zenith, counter, en))

dagman = pycondor.Dagman('tauRunner_cosmogenic_flux', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
