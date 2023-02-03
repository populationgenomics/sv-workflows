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
def main(): 
    b = hb.Batch("name", default_python_image=config['workflow']['driver_image'])
    j = b.new_python_job(name = "potato")

    text = (j.call(hello_world, 'alice')).as_str()
    #result = j.call(upper, hello_str)
    b.write_output(text, output_path('output/hello-alice.txt'))
    b.run(wait = False)
if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter