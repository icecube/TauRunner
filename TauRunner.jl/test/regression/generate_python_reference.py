#!/usr/bin/env python3
"""
Generate reference data from Python TauRunner for regression testing.

This script runs the Python TauRunner at various energies and angles,
saving results for comparison with the Julia implementation.

Usage:
    python generate_python_reference.py --output reference_data.json
"""

import argparse
import json
import sys
import os

# Add taurunner to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import numpy as np
from taurunner.body.earth import construct_earth
from taurunner.track.chord import Chord
from taurunner.track.radial import Radial
from taurunner.cross_sections import CrossSections
from taurunner.main import run_MC
from taurunner.utils import units


def generate_geometry_tests(earth):
    """Generate reference data for geometry calculations."""
    results = {
        'earth_radius_km': float(earth.radius / units.km),
        'density_profile': [],
        'chord_tests': [],
        'column_depth_tests': [],
    }

    # Density at various radii
    for r in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        rho = earth.get_density(r) / (units.gr / units.cm**3)
        results['density_profile'].append({
            'r': float(r),
            'density_gcm3': float(rho)
        })

    # Chord properties and column depths at various angles
    for theta_deg in [0, 10, 20, 30, 45, 60, 70, 80, 85, 89]:
        theta_rad = np.radians(theta_deg)
        chord = Chord(theta=theta_rad, depth=0.0)

        # Track length
        total_length = chord._l1 + chord._l2  # in units of radius

        # Column depth
        column_depth = chord.total_column_depth(earth)
        column_depth_gcm2 = column_depth / (units.gr / units.cm**2)

        # x_to_r at various positions
        x_to_r_values = []
        for x in [0.0, 0.25, 0.5, 0.75, 1.0]:
            r = chord.x_to_r(x)
            x_to_r_values.append({'x': float(x), 'r': float(r)})

        results['chord_tests'].append({
            'theta_deg': int(theta_deg),
            'total_length': float(total_length),
            'column_depth_gcm2': float(column_depth_gcm2),
            'x_to_r': x_to_r_values
        })

    return results


def generate_mc_tests(earth, xs, seed=42):
    """Generate reference data for Monte Carlo propagation."""
    results = []

    # Test configurations: (energy_eV, theta_deg, n_events)
    test_configs = [
        (1e15, 85, 100),   # 1 PeV, grazing
        (1e15, 60, 100),   # 1 PeV, moderate
        (1e15, 30, 100),   # 1 PeV, steep
        (1e17, 85, 100),   # 100 PeV, grazing
        (1e17, 60, 100),   # 100 PeV, moderate
        (1e18, 89, 100),   # 1 EeV, very grazing
        (1e18, 85, 100),   # 1 EeV, grazing
    ]

    for energy, theta_deg, n_events in test_configs:
        np.random.seed(seed)

        energies = np.ones(n_events) * energy
        thetas = np.ones(n_events) * np.radians(theta_deg)

        output = run_MC(
            energies,
            thetas,
            earth,
            xs,
            seed=seed,
            no_secondaries=False
        )

        # Extract tau neutrino results
        nutau_mask = np.abs(output['PDG_Encoding']) == 16
        survived_mask = output['Eout'] > 0

        tau_mask = np.abs(output['PDG_Encoding']) == 15

        result = {
            'energy_eV': energy,
            'theta_deg': theta_deg,
            'n_events': n_events,
            'seed': seed,
            # Nu_tau statistics
            'nutau_survived': int(np.sum(nutau_mask & survived_mask)),
            'nutau_mean_energy_out': float(np.mean(output['Eout'][nutau_mask & survived_mask])) if np.any(nutau_mask & survived_mask) else 0,
            'nutau_median_energy_out': float(np.median(output['Eout'][nutau_mask & survived_mask])) if np.any(nutau_mask & survived_mask) else 0,
            'nutau_std_energy_out': float(np.std(output['Eout'][nutau_mask & survived_mask])) if np.any(nutau_mask & survived_mask) else 0,
            # Tau statistics
            'tau_exited': int(np.sum(tau_mask & survived_mask)),
            'tau_mean_energy_out': float(np.mean(output['Eout'][tau_mask & survived_mask])) if np.any(tau_mask & survived_mask) else 0,
            # Interaction counts
            'mean_n_cc': float(np.mean(output['nCC'])),
            'mean_n_nc': float(np.mean(output['nNC'])),
            # Per-event data for exact comparison (first 10 events)
            'per_event': [
                {
                    'Eini': float(output['Eini'][i]),
                    'Eout': float(output['Eout'][i]),
                    'PDG': int(output['PDG_Encoding'][i]),
                    'nCC': int(output['nCC'][i]),
                    'nNC': int(output['nNC'][i]),
                }
                for i in range(min(10, n_events))
            ]
        }

        results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(description='Generate Python TauRunner reference data')
    parser.add_argument('--output', '-o', default='reference_data.json',
                        help='Output JSON file')
    parser.add_argument('--skip-mc', action='store_true',
                        help='Skip Monte Carlo tests (geometry only)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for MC')
    args = parser.parse_args()

    print("Constructing Earth model...")
    earth = construct_earth()

    print("Generating geometry reference data...")
    geometry_data = generate_geometry_tests(earth)

    reference = {
        'version': 'python',
        'geometry': geometry_data,
        'monte_carlo': []
    }

    if not args.skip_mc:
        print("Loading cross-sections...")
        xs = CrossSections('CSMS')

        print("Running Monte Carlo tests...")
        mc_data = generate_mc_tests(earth, xs, seed=args.seed)
        reference['monte_carlo'] = mc_data

    print(f"Writing reference data to {args.output}...")
    with open(args.output, 'w') as f:
        json.dump(reference, f, indent=2)

    print("Done!")


if __name__ == '__main__':
    main()
