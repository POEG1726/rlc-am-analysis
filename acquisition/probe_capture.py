
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import pyvisa, time
SCOPE=_bench.scope_address; GEN=_bench.generator_address
rm=pyvisa.ResourceManager(); sc=rm.open_resource(SCOPE); gn=rm.open_resource(GEN)
sc.timeout=120000; gn.timeout=5000; sc.chunk_size=4*1024*1024
try:
    # User-authorized diagnostic change: pure sine on generator CH2 only.
    gn.write('C2:MDWV STATE,OFF'); time.sleep(.2)
    st=gn.query('C2:MDWV?').strip(); print('MDWV',st)
    if 'STATE,OFF' not in st.upper(): raise RuntimeError('AM did not turn off')
    gn.write('C2:BSWV WVTP,SINE,FRQ,495000HZ,AMP,5V,OFST,0V')
    print('BSWV',gn.query('C2:BSWV?').strip())
    # Short diagnostic acquisition: 2 ms at 50 MSa/s, 100 kpts.
    sc.write(':TIM:SCAL 2E-4'); sc.write(':ACQ:MDEP 100K'); time.sleep(.5)
    for q in (':ACQ:POIN?',':ACQ:SRAT?',':ACQ:MDEP?',':TIM:SCAL?',':TRIG:STAT?'):
        print(q,sc.query(q).strip())
    for ch in (1,2): sc.write(f':WAV:SOUR C{ch}'); print('PRE',ch,sc.query_ascii_values(':WAV:PRE?') if False else 'binary descriptor below')
    sc.write(':TRIG:RUN'); time.sleep(.3); sc.write(':TRIG:STOP')
    for ch in (1,2):
        sc.write(f':WAV:SOUR C{ch}'); sc.write(':WAV:WIDT WORD'); sc.write(':WAV:BYT LSB'); sc.write(':WAV:INT 1'); sc.write(':WAV:STAR 0'); sc.write(':WAV:POIN 100000')
        sc.write(':WAV:PRE?'); b=sc.read_raw(); p=b.find(b'#'); n=int(chr(b[p+1])); k=int(b[p+2:p+2+n]); desc=b[p+2+n:p+2+n+k]
        import struct
        print('CH',ch,'desc',len(desc),'vdiv',struct.unpack_from('<f',desc,156)[0],'offset',struct.unpack_from('<f',desc,160)[0],'codes/div',struct.unpack_from('<f',desc,164)[0],'dt',struct.unpack_from('<f',desc,176)[0])
        sc.write(':WAV:DATA?'); b=sc.read_raw(); p=b.find(b'#'); n=int(chr(b[p+1])); k=int(b[p+2:p+2+n]); data=b[p+2+n:p+2+n+k]; print('DATA',len(data),'pts',len(data)//2)
    sc.write(':TRIG:RUN')
finally:
    sc.close();gn.close();rm.close()
