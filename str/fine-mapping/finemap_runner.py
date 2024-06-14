#!/usr/bin/env python3
"""
analysis-runner --dataset "bioheart" --access-level "test" --description "Run finemap" --output-dir "str/finemap" finemap_runner.py

"""
from cpg_utils.hail_batch import get_batch

def main():
    b = get_batch()
    input_file_path = 'gs://cpg-bioheart-test/str/associatr/finemap/example/data'
    data_in = b.read_input_group(**{
            'k':input_file_path+ '.k',
            'ld':input_file_path+ '.ld',
            'z' : input_file_path+ '.z',
            'master' : input_file_path+ '.master'})



    eh_job = b.new_job(name=f'Run finemap')
    eh_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/finemap:1.4.2')
    eh_job.command(
                    f"""
                # Create a temporary file
                temp_file=$(mktemp)


                # Write the required format to the file
                echo 'z;ld;snp;config;cred;log;k;n_samples' > $temp_file
                echo "{data_in.z};{data_in.ld};data.snp;data.config;data.cred;data.log;{data_in.k};5363" >> $temp_file
                finemap --sss --in-files temp_file --dataset 1
                """,
                )
    b.run(wait=False)
if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter