#!/usr/bin/env python3

"""

This script will run SusieR, fine-mapping tool using the residualized Y and X files (direct genotypes).

Required inputs:
- output from`files_prep_residualizer.py` (residualized Y and X files) (which requires files_prep_dosages.py to run first).
- List of eGenes to run SusieR on (ie Table S1 from the manuscript).

analysis-runner --dataset tenk10k --access-level test --memory 8G  --image "australia-southeast1-docker.pkg.dev/cpg-common/images/r-meta:r-fast" --description "Run SusiE MCV" --output-dir str/associatr/final_freeze/fine_mapping/susie_mcv/output_files \
susie_mcv_runner.py --table-s1-path=gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv --residualized-dir=gs://cpg-tenk10k-test/str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files/residualized  --max-parallel-jobs=1000 \
--chromosomes=chr9,chr10,chr13,chr14,chr15,chr16 --susie-cpu=0.5
"""

import click
import math
import hailtop.batch as hb
import pandas as pd
import gc

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def susie_runner(input_dir, df_batch, cell_type, chrom, num_causal_variants, num_iterations):
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(susieR)')
    ro.r('library(tidyverse)')

    for gene in df_batch['gene_name'].unique():
        if to_path(
                output_path(f'{cell_type}/{chrom}/{gene}_100kb.tsv', 'analysis'),
            ).exists():
            print(f'SusieR file for {gene} in {cell_type} already exists, skipping.')
            continue
        print(f'Running SuSiE for {gene} in {cell_type} on {chrom}')

        # load in X and y file paths
        X_df = pd.read_csv(f'{input_dir}/{cell_type}/{chrom}/{gene}_{cell_type}_meta_cleaned_X_resid.csv')
        y_df = pd.read_csv(f'{input_dir}/{cell_type}/{chrom}/{gene}_{cell_type}_meta_cleaned_y_resid.csv')

        with (ro.default_converter + pandas2ri.converter).context():
            X = ro.conversion.get_conversion().py2rpy(X_df)

        with (ro.default_converter + pandas2ri.converter).context():
            y = ro.conversion.get_conversion().py2rpy(y_df)

        # === 4. Assign to R environment ===
        ro.globalenv['X'] = X
        ro.globalenv['y'] = y
        ro.globalenv['num_causal_variants'] = num_causal_variants
        ro.globalenv['num_iterations'] = num_iterations

        # === 5. Run SuSiE ===
        ro.r(
            '''
        X <- subset(X, select = -sample)
        x_input = as.matrix(X)
        variant_ids <- colnames(x_input)
        y = subset(y, select=-sample)
        y_input = y[,1]
        # run susieR
        susie_fit <- susie(x_input, y_input, L = num_causal_variants, max_iter = num_iterations)

        #capture susie output results to save later
        raw_output = capture.output(summary(susie_fit))
            '''
        )

        # convert raw output to python
        raw_output_python = ro.r('raw_output')

        # write raw output to GCS
        with to_path(output_path(f"{cell_type}/{chrom}/{gene}_100kb_output.txt", 'analysis')).open('w') as file:
            file.write(str(raw_output_python))

        ro.r('coord_df <- data.frame(variant_id = variant_ids)')
        ro.r('coord_df$chr <- sub("\\\\..*", "", coord_df$variant_id)')
        print('extracted coordinates from variant IDs')

        ro.r(
            '''
        # Get credible sets from SuSiE fit
        susie_cs <- susie_get_cs(susie_fit, X = x_input)

        # Initialize annotation vectors
        n_variants <- length(variant_ids)
        cs_id <- rep(NA_real_, n_variants)
        cs_size <- rep(NA_real_, n_variants)
        max_pip_in_cs <- rep(NA_real_, n_variants) '''
        )
        print('obtained credible sets and initialized vectors')

        ro.r(
            '''
        # Annotate variants with CS membership
        for (i in seq_along(susie_cs$cs)) {
        idx <- susie_cs$cs[[i]]
        cs_id[idx] <- i
        cs_size[idx] <- length(idx)
        max_pip_in_cs[idx] <- max(susie_fit$pip[idx], na.rm = TRUE)
        }
        '''
        )
        print('annotated variants with CS membership')

        ro.r(
            '''
        # Combine results into final dataframe
        final_df <- data.frame(
        variant_id = variant_ids,
        pip = susie_fit$pip,
        pip_prune = susie_get_pip(susie_fit, prune_by_cs = TRUE),
        cs_id = cs_id,
        cs_size = cs_size,
        max_pip_in_cs = max_pip_in_cs
        ) |>
        merge(coord_df, by = "variant_id") |>
        dplyr::arrange(pip)

        '''
        )
        print('combined results into final dataframe')

        with (ro.default_converter + pandas2ri.converter).context():
            susie_output_df = ro.conversion.get_conversion().rpy2py(ro.r('final_df'))
        print('converted back to pandas df')

        # write dataframe to GCS
        susie_output_df.to_csv(
            output_path(f"{cell_type}/{chrom}/{gene}_100kb.tsv", 'analysis'),
            sep='\t',
            index=False,
        )
        del X_df, y_df, X, y,susie_output_df
        gc.collect()


@click.option('--table-s1-path', help='Table S1 with eGenes to run SusieR on')
@click.option('--residualized-dir', help='Directory with residualized Y and X files')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs', default=500)
@click.option('--num_iterations', help='Number of iterations for SusieR', default=100)
@click.option('--susie-cpu', help='CPU for SusieR job', default=0.25)
@click.option('--num-causal-variants', help='Number of causal variants to estimate', default=10)
@click.option('--always-run', help='Job set to always run', is_flag=True)
@click.option('--chromosomes', help='Chromosomes to process, comma separated', default='chr21')
@click.command()
def main(
    table_s1_path,
    residualized_dir,
    max_parallel_jobs,
    num_iterations,
    susie_cpu,
    num_causal_variants,
    always_run,
    chromosomes
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
    df = df.drop_duplicates(subset=['cell_type', 'gene_name', 'chr'])
    for chrom in chromosomes.split(','):
        df_chrom = df[df['chr'] == chrom]
        for cell_type in df_chrom['cell_type'].unique():
            df_cell = df_chrom[df_chrom['cell_type'] == cell_type]
            # Split df_cell into batches of 30 rows
            num_batches = math.ceil(len(df_cell) / 30)
            for i in range(num_batches):
                df_batch = df_cell.iloc[i*30 : (i+1)*30]
                susie_job = b.new_python_job(
                f'SusieR {cell_type} {chrom} batch {i}',
                )
                susie_job.cpu(susie_cpu)
                if always_run:
                    susie_job.always_run()
                susie_job.call(
                    susie_runner,
                    residualized_dir,
                    df_batch,
                    cell_type,
                    chrom,
                    num_causal_variants,
                    num_iterations,
                )
                manage_concurrency_for_job(susie_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()
