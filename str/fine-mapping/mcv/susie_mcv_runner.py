#!/usr/bin/env python3

"""

This script will run SusieR, fine-mapping tool using the residualized Y and X files (direct genotypes).

Required inputs:
- output from`files_prep_residualizer.py` (residualized Y and X files) (which requires files_prep_dosages.py to run first).
- List of eGenes to run SusieR on (ie Table S1 from the manuscript).

analysis-runner --dataset bioheart --access-level test --memory 4G  --image "australia-southeast1-docker.pkg.dev/cpg-common/images/r-meta:susie" --description "Run SusiE MCV" --output-dir tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/output_files \
susie_mcv_runner.py --table-s1-path=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/TableS1.csv --residualized-dir=gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files
"""

import click

import hailtop.batch as hb
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def susie_runner(input_dir, gene, cell_type, num_causal_variants, num_iterations):
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(susieR)')
    ro.r('library(tidyverse)')

    # load in X and y file paths
    X_df = pd.read_csv(f'{input_dir}/{gene}_{cell_type}_meta_cleaned_X_resid.csv')
    y_df = pd.read_csv(f'{input_dir}/{gene}_{cell_type}_meta_cleaned_y_resid.csv')

    with (ro.default_converter + pandas2ri.converter).context():
        X = ro.conversion.get_conversion().py2rpy(X_df)

    with (ro.default_converter + pandas2ri.converter).context():
        y = ro.conversion.get_conversion().py2rpy(y_df)

    # === 4. Assign to R environment ===
    ro.globalenv['X'] = X
    ro.globalenv['y'] = y

    # === 5. Run SuSiE ===
    ro.r('''
    library(susieR)
    X <- subset(X, select = -sample)
    x_input = as.matrix(X)
    variant_ids <- colnames(x_input)
    y = subset(y, select=-sample)
    y_input = y[,1]
    # run susieR
    susie_fit <- susie(x_input, y_input, L = 10)

    #capture susie output results to save later
    raw_output = capture.output(summary(susie_fit))

    coord_df <- data.frame(variant_id = variant_ids)
    coord_df$chr <- sub("\\..*", "", coord_df$variant_id)
    # Extract position (the number between the first and second dot)
    coord_df$pos <- as.integer(sub("^[^\\.]+\\.([^\\.]+)\\..*$", "\\1", coord_df$variant_id))

    susie_cs <- susie_get_cs(susie_fit, X = x_input)

    # Initialize annotation vectors
    n_variants <- length(variant_ids)
    cs_id <- rep(NA_integer_, n_variants)
    cs_size <- rep(NA_integer_, n_variants)
    max_pip_in_cs <- rep(NA_real_, n_variants)

    # Loop over credible sets
    for (i in seq_along(susie_cs$cs)) {
    cs_variants <- susie_cs$cs[[i]]
    purity <- susie_cs$purity[[i]]
    max_pip <- max(susie_fit$pip[cs_variants], na.rm = TRUE)

    cs_id[cs_variants] <- i
    cs_size[cs_variants] <- length(cs_variants)
    max_pip_in_cs[cs_variants] <- max_pip
    }


    # -----------------------------
    # Step 5: Compile final output
    # -----------------------------
    pip_df <- data.frame(
    variant_id = variant_ids,
    pip = susie_fit$pip,
    pip_prune = susie_get_pip(susie_fit, prune_by_cs = TRUE),
    cs_id = cs_id,
    cs_size = cs_size,
    max_pip_in_cs = max_pip_in_cs
    )

    final_df <- merge(pip_df, coord_df, by = "variant_id")
    final_df <- final_df[order(final_df$chr, final_df$pos), ]




    ''')

    # convert raw output to python
    raw_output_python = ro.r('raw_output')

    # write raw output to GCS
    with to_path(output_path(f"{cell_type}/{gene}_100kb_output.txt", 'analysis')).open('w') as file:
        file.write(str(raw_output_python))

    with (ro.default_converter + pandas2ri.converter).context():
        susie_output_df = ro.conversion.get_conversion().rpy2py(ro.r('final_df'))
    print('converted back to pandas df')

    # write dataframe to GCS
    susie_output_df.to_csv(
        output_path(f"susie/{cell_type}/{gene}_100kb.tsv", 'analysis'),
        sep='\t',
        index=False,
    )


@click.option('--table-s1-path', help='Table S1 with eGenes to run SusieR on')
@click.option('--residualized-dir', help='Directory with residualized Y and X files')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs', default=500)
@click.option('--num_iterations', help='Number of iterations for SusieR', default=100)
@click.option('--susie-cpu', help='CPU for SusieR job', default=8)
@click.option('--num-causal-variants', help='Number of causal variants to estimate', default=10)
@click.option('--always-run', help='Job set to always run', is_flag=True)
@click.command()
def main(
    table_s1_path,
    residualized_dir,
    max_parallel_jobs,
    num_iterations,
    susie_cpu,
    num_causal_variants,
    always_run,
):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name='Run susieR')

    df = pd.read_csv(table_s1_path)
    df = df[df['chr'] == 'chr1']  # For testing, only run on chr1
    df = df.drop_duplicates(subset=['cell_type', 'gene_name'])
    for row in df.itertuples():
        cell_type = row.cell_type
        gene = row.gene_name
        if to_path(
            output_path(f'susie/{cell_type}/{gene}_100kb.tsv', 'analysis'),
        ).exists():
            print(f'SusieR file for {gene} in {cell_type} already exists, skipping.')
            continue
        susie_job = b.new_python_job(
            f'SusieR for {gene}: {cell_type}',
        )
        susie_job.cpu(susie_cpu)
        if always_run:
            susie_job.always_run()
        susie_job.call(
            susie_runner,
            residualized_dir,
            gene,
            cell_type,
            num_causal_variants,
            num_iterations,
        )
        break  # For testing, only run one gene at a time
        manage_concurrency_for_job(susie_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()
