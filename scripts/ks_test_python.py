#!/usr/bin/env python3
"""Run tau neutrinos through Earth at several angles and save outgoing energies."""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from taurunner.body import construct_earth
from taurunner.cross_sections import CrossSections
from taurunner.main import run_MC

def main():
    n_events = int(sys.argv[1])
    thetas_deg = [float(x) for x in sys.argv[2].split(",")]
    energy = float(sys.argv[3])
    outdir = sys.argv[4]
    seed = int(sys.argv[5])

    earth = construct_earth()
    xs = CrossSections("CSMS")

    for theta_deg in thetas_deg:
        theta_rad = np.deg2rad(theta_deg)
        energies = np.full(n_events, energy)
        thetas = np.full(n_events, theta_rad)

        results = run_MC(energies, thetas, earth, xs,
                         seed=seed, flavor=16,
                         no_secondaries=True, no_losses=False)

        # Collect final energies for particles that exited (primaries only since no_secondaries=True)
        exited = results["final_position"] >= 1.0
        out_energies = results["Eout"][exited]

        fname = os.path.join(outdir, f"python_theta_{theta_deg}_n_{n_events}.csv")
        np.savetxt(fname, out_energies, delimiter=",")
        print(f"θ={theta_deg}°: {len(out_energies)}/{n_events} exited, saved to {fname}")

if __name__ == "__main__":
    main()
