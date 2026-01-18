#!/usr/bin/env python3
"""
Convert TauRunner Python pickle files to Julia-compatible JLD2/JSON format.

This script extracts cross-section spline data from Python pickle files
and saves them in a format that Julia can read.

Usage:
    python convert_pickle_to_jld2.py --input /path/to/taurunner/resources --output /path/to/TauRunner.jl/data
"""

import argparse
import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np


def extract_spline_data(spline):
    """Extract data from a scipy spline object."""
    if hasattr(spline, 'get_knots'):
        # UnivariateSpline
        return {
            'type': '1d',
            'knots': spline.get_knots().tolist(),
            'coefficients': spline.get_coeffs().tolist(),
            'degree': int(spline.k) if hasattr(spline, 'k') else 3,
        }
    elif hasattr(spline, 'x') and hasattr(spline, 'y'):
        # interp1d or similar
        return {
            'type': '1d',
            'knots': np.array(spline.x).tolist(),
            'values': np.array(spline.y).tolist(),
        }
    elif isinstance(spline, np.ndarray):
        # Raw numpy array
        return {
            'type': 'array',
            'data': spline.tolist(),
        }
    else:
        # Try to extract what we can
        return {
            'type': 'unknown',
            'repr': repr(spline)[:200],
        }


def convert_cross_section_pickle(input_path, output_path):
    """Convert a single cross-section pickle file."""
    try:
        with open(input_path, 'rb') as f:
            data = pickle.load(f)

        extracted = extract_spline_data(data)

        # Save as JSON for Julia to read
        with open(output_path, 'w') as f:
            json.dump(extracted, f, indent=2)

        return True
    except Exception as e:
        print(f"Error converting {input_path}: {e}", file=sys.stderr)
        return False


def convert_npy_file(input_path, output_path):
    """Convert a numpy .npy file to JSON."""
    try:
        data = np.load(input_path, allow_pickle=True)

        if isinstance(data, np.ndarray):
            output = {
                'type': 'array',
                'shape': list(data.shape),
                'data': data.tolist(),
            }
        else:
            output = {
                'type': 'unknown',
                'data': str(data)[:200],
            }

        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)

        return True
    except Exception as e:
        print(f"Error converting {input_path}: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Convert TauRunner Python data files to Julia format'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Path to taurunner/resources directory'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Path to TauRunner.jl/data directory'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )

    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)

    if not input_dir.exists():
        print(f"Input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(1)

    # Create output directory structure
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert cross-section tables
    xs_input = input_dir / 'cross_section_tables'
    if xs_input.exists():
        # Files are in the cross_section_tables directory with prefixes
        for model in ['dipole', 'CSMS']:
            model_output = output_dir / 'cross_sections' / model
            model_output.mkdir(parents=True, exist_ok=True)

            # Look for files with model prefix (e.g., CSMS_nu_p_sigma_CC.pkl)
            for pkl_file in xs_input.glob(f'{model}_*.pkl'):
                # Remove model prefix from output filename
                out_name = pkl_file.stem.replace(f'{model}_', '')
                out_file = model_output / (out_name + '.json')
                if args.verbose:
                    print(f"Converting {pkl_file} -> {out_file}")
                convert_cross_section_pickle(pkl_file, out_file)

        # Also convert non-prefixed pkl files to a 'default' directory
        default_output = output_dir / 'cross_sections' / 'default'
        default_output.mkdir(parents=True, exist_ok=True)

        for pkl_file in xs_input.glob('*.pkl'):
            if not pkl_file.stem.startswith(('CSMS_', 'dipole_')):
                out_file = default_output / (pkl_file.stem + '.json')
                if args.verbose:
                    print(f"Converting {pkl_file} -> {out_file}")
                convert_cross_section_pickle(pkl_file, out_file)

        # Convert npy files
        for npy_file in xs_input.glob('*.npy'):
            out_file = default_output / (npy_file.stem + '.json')
            if args.verbose:
                print(f"Converting {npy_file} -> {out_file}")
            convert_npy_file(npy_file, out_file)

    # Convert secondaries splines
    sec_input = input_dir / 'secondaries_splines'
    if sec_input.exists():
        sec_output = output_dir / 'secondaries'
        sec_output.mkdir(parents=True, exist_ok=True)

        for pkl_file in sec_input.glob('*.pkl'):
            out_file = sec_output / (pkl_file.stem + '.json')
            if args.verbose:
                print(f"Converting {pkl_file} -> {out_file}")
            convert_cross_section_pickle(pkl_file, out_file)

        for npy_file in sec_input.glob('*.npy'):
            out_file = sec_output / (npy_file.stem + '.json')
            if args.verbose:
                print(f"Converting {npy_file} -> {out_file}")
            convert_npy_file(npy_file, out_file)

    # Copy solar models (text files can be used directly)
    solar_input = input_dir / 'solar_models'
    if solar_input.exists():
        solar_output = output_dir / 'solar_models'
        solar_output.mkdir(parents=True, exist_ok=True)

        for txt_file in solar_input.glob('*.txt'):
            import shutil
            if args.verbose:
                print(f"Copying {txt_file}")
            shutil.copy(txt_file, solar_output / txt_file.name)

    # Convert GZK flux file
    gzk_file = input_dir / 'gzk_cdf_phi_spline.npy'
    if gzk_file.exists():
        out_file = output_dir / 'gzk_cdf_phi_spline.json'
        if args.verbose:
            print(f"Converting {gzk_file} -> {out_file}")
        convert_npy_file(gzk_file, out_file)

    print("Conversion complete!")


if __name__ == '__main__':
    main()
