"""Configure and run one preserved bench script. Listing/checking never opens VISA."""
import argparse
import json
import os
from pathlib import Path
import runpy
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(p.name for p in Path(__file__).parent.glob('*.py')
                 if p.name not in {'run.py', 'acquisition_config.py'})


def configuration(path):
    config = json.loads(path.read_text(encoding='utf-8-sig'))
    for key in ('scope_address', 'generator_address'):
        if not isinstance(config.get(key), str) or not config[key].strip():
            raise ValueError(f'Set {key} in the private configuration file.')
    workspace = (ROOT / config.get('workspace', 'local_runs/bench01')).resolve()
    if not workspace.is_relative_to(ROOT / 'local_runs') or workspace == ROOT / 'local_runs':
        raise ValueError('workspace must be a subdirectory of local_runs/.')
    return config, workspace


def seed_reference(workspace, script):
    """Copy only explicit prerequisite inputs; never seed a resume-capable result table."""
    reference = ROOT / 'data/characterization'
    inputs = {
        'diagnose_am_1khz_low_drive.py': ['vertical_scale_corrections.json'],
        'sweep_low_drive_snr.py': ['vertical_scale_corrections.json'],
        'verify_low_drive_repeatability.py': ['vertical_scale_corrections.json'],
        'verify_low_drive_repeatability_bab.py': ['vertical_scale_corrections.json'],
        'refine_resonance.py': ['resonance_scan.csv','characterization_parameters.json'],
        'compare_drive_levels.py': ['characterization_parameters.json'],
        'summarize_drive_amplitude.py': ['scope_vertical_range_calibration.csv',
            'drive_amplitude_frequency_scan.csv','drive_amplitude_sweep.csv','drive_amplitude_harmonics.csv'],
    }
    for name in inputs.get(script, []):
        target = workspace / name
        if not target.exists():
            shutil.copyfile(reference / name, target)
    if script == 'collect_formal_am_full_rerun.py':
        target = workspace / 'formal_data/processed/formal_am_summary.csv'
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / 'data/am_initial/formal_am_summary.csv', target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('script', nargs='?', choices=SCRIPTS)
    parser.add_argument('--list', action='store_true', help='List scripts without connecting to equipment.')
    parser.add_argument('--check', action='store_true', help='Validate local config only; no equipment access.')
    parser.add_argument('--config', type=Path, default=ROOT / 'config.local.json')
    args = parser.parse_args()
    if args.list:
        print('\n'.join(SCRIPTS))
        return
    if not args.script and not args.check:
        parser.error('Choose a script or use --list / --check.')
    try:
        config, workspace = configuration(args.config)
    except (OSError, ValueError, TypeError) as exc:
        parser.error('Local configuration is missing or invalid. Copy config.example.json to '
                     'config.local.json and configure both instruments and workspace. ' +
                     (str(exc) if isinstance(exc, ValueError) else ''))
    if args.check:
        print('Configuration valid. No VISA session opened; resource addresses withheld.')
        return
    script = Path(__file__).parent / args.script
    # The extreme acquisition is a historical resume routine, not a new-run entry point.
    if args.script == 'collect_drive_dependence_extreme.py':
        required = workspace / 'drive_dependence_extreme_12Vpp_20260930_run02/drive_dependence_results.csv'
        if not required.exists():
            parser.error('This archived script resumes a private partial run; see docs/ACQUISITION.md.')
    workspace.mkdir(parents=True, exist_ok=True)
    seed_reference(workspace, args.script)
    os.environ.update(RLC_ACQUISITION_RUNNER='1', RLC_SCOPE_ADDRESS=config['scope_address'],
                      RLC_GENERATOR_ADDRESS=config['generator_address'], RLC_WORKSPACE=str(workspace))
    if config.get('visa_library'):
        import pyvisa
        original = pyvisa.ResourceManager
        library = config['visa_library']
        pyvisa.ResourceManager = lambda *a, **kw: original(*a, **kw) if a or kw else original(library)
    sys.path.insert(0, str(script.parent))
    os.chdir(workspace)  # Also redirects relative filenames in older scripts.
    runpy.run_path(str(script), run_name='__main__')


if __name__ == '__main__':
    main()
