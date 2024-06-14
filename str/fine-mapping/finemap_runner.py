#!/usr/bin/env python3
"""
analysis-runner --dataset "bioheart" --access-level "test" --description "Run finemap" --output-dir "str/finemap" finemap_runner.py

"""
from cpg_utils.hail_batch import get_batch, image_path, output_path, reference_path

def main():
    b = get_batch()
    eh_job = b.new_job(name=f'Run finemap')
    eh_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/finemap:1.4.2')
    eh_job.command(
                    f"""
                finemap
                """,
                )
    b.run(wait=False)
if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter