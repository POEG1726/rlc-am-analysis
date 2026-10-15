
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import pyvisa

SCOPE_ADDR = _bench.scope_address
GEN_ADDR   = _bench.generator_address

rm = pyvisa.ResourceManager()

scope = rm.open_resource(SCOPE_ADDR)
gen = rm.open_resource(GEN_ADDR)

scope.timeout = 5000
gen.timeout = 5000

print("=== IDENTIFICATION ===")
print("Scope:    ", scope.query("*IDN?").strip())
print("Generator:", gen.query("*IDN?").strip())

print("\n=== GENERATOR CH2 ===")

print("\nBasic waveform:")
print(gen.query("C2:BSWV?").strip())

print("\nModulation:")
print(gen.query("C2:MDWV?").strip())

print("\nOutput:")
print(gen.query("C2:OUTP?").strip())

scope.close()
gen.close()
rm.close()
