
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import pyvisa, json, time
SCOPE=_bench.scope_address
GEN=_bench.generator_address
rm=pyvisa.ResourceManager()
sc=rm.open_resource(SCOPE); gn=rm.open_resource(GEN)
sc.timeout=5000; gn.timeout=5000
sc.write_termination='\n'; sc.read_termination='\n'
gn.write_termination='\n'; gn.read_termination='\n'
print('SCOPE',sc.query('*IDN?').strip())
print('GEN',gn.query('*IDN?').strip())
queries=[':ACQ:POIN?', ':ACQ:SRAT?', ':ACQ:MDEP?', ':TIM:SCAL?', ':TRIG:STAT?', ':WAV:PRE?']
for q in queries:
    try:
        ans=sc.query(q).strip()
        print('SCOPE',q,ans[:240])
    except Exception as e: print('SCOPE_ERROR',q,repr(e))
for ch in (1,2):
    for q in (f':CHAN{ch}:SCAL?',f':CHAN{ch}:OFFS?',f':CHAN{ch}:PROB?',f':CHAN{ch}:COUP?',f':CHAN{ch}:IMP?',f':CHAN{ch}:BWL?'):
        try: print('SCOPE',q,sc.query(q).strip())
        except Exception as e: print('SCOPE_ERROR',q,repr(e))
print('GEN BSWV',gn.query('C2:BSWV?').strip())
print('GEN MDWV',gn.query('C2:MDWV?').strip())
print('GEN OUTP',gn.query('C2:OUTP?').strip())
sc.close();gn.close();rm.close()
