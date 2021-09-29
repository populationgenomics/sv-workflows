#!/usr/bin/env python

from sample_metadata.api import SampleApi
sapi = SampleApi()
s = sapi.get_sample_id_map_by_internal(["CPG3293", "CPG3301"])
print(s)
