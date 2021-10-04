#!/usr/bin/env python

from sample_metadata.api import SampleApi
sapi = SampleApi()
cpg_samples = ["CPG10017", "CPG10025"]
tob = sapi.get_sample_id_map_by_internal(cpg_samples)
print(f'Got {len(cpg_samples)} samples:')
print(tob)
