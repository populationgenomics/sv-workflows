#!/usr/bin/env python3
import os
import logging


from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch
from google.cloud import storage
config = get_config()

def hello_world(name):
    return f'hello {name}'


def upper(s):
    return s.upper()


b = get_batch()
j = b.new_python_job(name = "potato")
hello_str = j.call(hello_world, 'alice')
result = j.call(upper, hello_str)
b.write_output(result.as_str(), 'output/hello-alice.txt')
b.run()