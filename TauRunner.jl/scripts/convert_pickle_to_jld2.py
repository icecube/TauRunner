#!/usr/bin/env python3
"""
Convert TauRunner Python pickle files to Julia-compatible JSON format.

This script extracts cross-section spline data from Python pickle files
and saves them in a format that Julia can read with Interpolations.jl.

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
    spline_type = type(spline).__name__

    # Check for 2D spline FIRST (RectBivariateSpline also has get_knots)
    if spline_type == 'RectBivariateSpline' or (hasattr(spline, 'tck') and hasattr(spline, 'degrees')):
        # 2D RectBivariateSpline
        tx, ty, c = spline.tck
        kx, ky = spline.degrees

        return {
            'type': 'RectBivariateSpline',
            'tx': np.array(tx).tolist(),  # x knots
            'ty': np.array(ty).tolist(),  # y knots
            'coefficients': np.array(c).tolist(),
            'kx': int(kx),
            'ky': int(ky),
        }

    elif spline_type == 'UnivariateSpline' or (hasattr(spline, 'get_knots') and hasattr(spline, 'get_coeffs')):
        # 1D UnivariateSpline
        knots = spline.get_knots()
        coeffs = spline.get_coeffs()
        # Get degree from internal data
        k = 3  # default
        if hasattr(spline, '_data') and len(spline._data) > 5:
            k = int(spline._data[5])

        return {
            'type': 'UnivariateSpline',
            'knots': knots.tolist(),
            'coefficients': coeffs.tolist(),
            'degree': k,
        }

    elif hasattr(spline, 'x') and hasattr(spline, 'y'):
        # interp1d or similar
        return {
            'type': 'interp1d',
            'x': np.array(spline.x).tolist(),
            'y': np.array(spline.y).tolist(),
        }

    elif isinstance(spline, np.ndarray):
        # Raw numpy array
        return {
            'type': 'array',
            'shape': list(spline.shape),
            'data': spline.tolist(),
        }

    else:
        # Try to extract what we can
        print(f"  Warning: Unknown spline type {spline_type}", file=sys.stderr)
        return {
            'type': 'unknown',
            'python_type': spline_type,
            'repr': repr(spline)[:500],
        }


def convert_cross_section_pickle(input_path, output_path):
    """Convert a single cross-section pickle file."""
    try:
        with open(input_path, 'rb') as f:
            try:
                data = pickle.load(f)
            except UnicodeDecodeError:
                # Python 2 pickle
                f.seek(0)
                data = pickle.load(f, encoding='latin1')

        extracted = extract_spline_data(data)

        # Save as JSON for Julia to read
        with open(output_path, 'w') as f:
            json.dump(extracted, f)

        return True
    except Exception as e:
        print(f"Error converting {input_path}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False


def convert_npy_file(input_path, output_path):
    """Convert a numpy .npy file to JSON."""
    try:
        # Try loading with different encodings for Python 2/3 compatibility
        try:
            data = np.load(input_path, allow_pickle=True)
        except (UnicodeDecodeError, ValueError):
            # Python 2 pickle - need fix_imports and latin1 encoding
            data = np.load(input_path, allow_pickle=True, fix_imports=True, encoding='latin1')

        # Check if it's a spline object stored in npy (0-d array with object)
        if isinstance(data, np.ndarray) and data.ndim == 0 and data.dtype == object:
            item = data.item()
            if hasattr(item, 'get_knots') or hasattr(item, 'tck'):
                output = extract_spline_data(item)
            else:
                output = {
                    'type': 'object',
                    'repr': repr(item)[:500],
                }
        elif isinstance(data, np.ndarray):
            output = {
                'type': 'array',
                'shape': list(data.shape),
                'data': data.tolist(),
            }
        else:
            output = {
                'type': 'unknown',
                'repr': str(data)[:500],
            }

        with open(output_path, 'w') as f:
            json.dump(output, f)

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

    success_count = 0
    fail_count = 0

    # Convert cross-section tables
    xs_input = input_dir / 'cross_section_tables'
    if xs_input.exists():
        # Files are in the cross_section_tables directory with prefixes
        for model in ['dipole', 'CSMS']:
            model_output = output_dir / 'cross_sections' / model.lower()
            model_output.mkdir(parents=True, exist_ok=True)

            # Look for files with model prefix (e.g., CSMS_nu_p_sigma_CC.pkl)
            for pkl_file in xs_input.glob(f'{model}_*.pkl'):
                # Remove model prefix from output filename
                out_name = pkl_file.stem.replace(f'{model}_', '')
                out_file = model_output / (out_name + '.json')
                if args.verbose:
                    print(f"Converting {pkl_file.name} -> {out_file}")
                if convert_cross_section_pickle(pkl_file, out_file):
                    success_count += 1
                else:
                    fail_count += 1

        # Also convert non-prefixed pkl files to a 'default' directory
        default_output = output_dir / 'cross_sections' / 'default'
        default_output.mkdir(parents=True, exist_ok=True)

        for pkl_file in xs_input.glob('*.pkl'):
            if not pkl_file.stem.startswith(('CSMS_', 'dipole_')):
                out_file = default_output / (pkl_file.stem + '.json')
                if args.verbose:
                    print(f"Converting {pkl_file.name} -> {out_file}")
                if convert_cross_section_pickle(pkl_file, out_file):
                    success_count += 1
                else:
                    fail_count += 1

        # Convert npy files
        for npy_file in xs_input.glob('*.npy'):
            out_file = default_output / (npy_file.stem + '.json')
            if args.verbose:
                print(f"Converting {npy_file.name} -> {out_file}")
            if convert_npy_file(npy_file, out_file):
                success_count += 1
            else:
                fail_count += 1

    # Convert secondaries splines
    sec_input = input_dir / 'secondaries_splines'
    if sec_input.exists():
        sec_output = output_dir / 'secondaries'
        sec_output.mkdir(parents=True, exist_ok=True)

        for pkl_file in sec_input.glob('*.pkl'):
            out_file = sec_output / (pkl_file.stem + '.json')
            if args.verbose:
                print(f"Converting {pkl_file.name} -> {out_file}")
            if convert_cross_section_pickle(pkl_file, out_file):
                success_count += 1
            else:
                fail_count += 1

        for npy_file in sec_input.glob('*.npy'):
            out_file = sec_output / (npy_file.stem + '.json')
            if args.verbose:
                print(f"Converting {npy_file.name} -> {out_file}")
            if convert_npy_file(npy_file, out_file):
                success_count += 1
            else:
                fail_count += 1

    # Convert tau decay splines
    decay_input = input_dir / 'tau_decay_tables'
    if decay_input.exists():
        decay_output = output_dir / 'tau_decay'
        decay_output.mkdir(parents=True, exist_ok=True)

        for pkl_file in decay_input.glob('*.pkl'):
            out_file = decay_output / (pkl_file.stem + '.json')
            if args.verbose:
                print(f"Converting {pkl_file.name} -> {out_file}")
            if convert_cross_section_pickle(pkl_file, out_file):
                success_count += 1
            else:
                fail_count += 1

        for npy_file in decay_input.glob('*.npy'):
            out_file = decay_output / (npy_file.stem + '.json')
            if args.verbose:
                print(f"Converting {npy_file.name} -> {out_file}")
            if convert_npy_file(npy_file, out_file):
                success_count += 1
            else:
                fail_count += 1

    # Copy solar models (text files can be used directly)
    solar_input = input_dir / 'solar_models'
    if solar_input.exists():
        solar_output = output_dir / 'solar_models'
        solar_output.mkdir(parents=True, exist_ok=True)

        import shutil
        for txt_file in solar_input.glob('*.txt'):
            if args.verbose:
                print(f"Copying {txt_file.name}")
            shutil.copy(txt_file, solar_output / txt_file.name)
            success_count += 1

    # Convert GZK flux file
    gzk_file = input_dir / 'gzk_cdf_phi_spline.npy'
    if gzk_file.exists():
        out_file = output_dir / 'gzk_cdf_phi_spline.json'
        if args.verbose:
            print(f"Converting {gzk_file.name} -> {out_file}")
        if convert_npy_file(gzk_file, out_file):
            success_count += 1
        else:
            fail_count += 1

    print(f"\nConversion complete! Success: {success_count}, Failed: {fail_count}")


if __name__ == '__main__':
    main()
