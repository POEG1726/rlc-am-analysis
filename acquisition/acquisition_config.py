"""Private runtime context for the archived bench scripts."""
from dataclasses import dataclass
import os
from pathlib import Path


@dataclass(frozen=True)
class BenchContext:
    scope_address: str
    generator_address: str
    workspace: Path


def require_context():
    if os.environ.get('RLC_ACQUISITION_RUNNER') != '1':
        raise RuntimeError('Use acquisition/run.py with a local configuration. See docs/ACQUISITION.md.')
    return BenchContext(os.environ['RLC_SCOPE_ADDRESS'], os.environ['RLC_GENERATOR_ADDRESS'],
                        Path(os.environ['RLC_WORKSPACE']))
