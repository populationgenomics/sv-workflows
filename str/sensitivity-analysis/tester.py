#!/usr/bin/env python3
import os
import logging


from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path
import hailtop.batch as hb


config = get_config()

def hello_world(name):
    return f'hello {name}'


def upper(s):
    return s.upper()
"""
backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
"""

b = hb.Batch("name", default_python_image=config['workflow']['driver_image'])
j = b.new_python_job(name = "potato")
j.declare_resource_group(
            eh_output={
                'txt': '{root}.txt'
            }
        )
j.eh_output = (j.call(hello_world, 'alice')).as_str()
#result = j.call(upper, hello_str)
b.write_output(j.eh_output, output_path('output/hello-alice.txt'))
b.run()