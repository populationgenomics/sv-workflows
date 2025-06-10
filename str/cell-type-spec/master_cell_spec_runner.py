#!/usr/bin/env python3

"""
This script performs cell type specificity analysis as per Cuomo et al. 2025 procedure.

analysis-runner --dataset tenk10k --access-level test --description "Cell type spec" --output-dir str/cell-spec-cuomo/cell_spec_results/v2 master_cell_spec_runner.py \
--meta-scen2-path=gs://cpg-tenk10k-test-analysis/str/cell-spec-cuomo/ab_effect_sizes/meta_results
"""

import numpy as np
import pandas as pd
import click
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def process_scenario_1_2(gene, chrom, cell_types_without_egene, associatr_path):
    """

    Processes scenario 1 and 2
    """
    scenario_1_2_dict = {}
    for cell_type in cell_types_without_egene:
        try:
            df = pd.read_csv(f'{associatr_path}/{cell_type}/{chrom}/{gene}_100000bp_meta_results.tsv', sep='\t')
            scenario_1_2_dict[cell_type] = 2
        except:
            scenario_1_2_dict[cell_type] = 1
    return scenario_1_2_dict


def process_same_tr(cell_types_with_same_eTR, row_cell_type, row_coeff):
    """

    Processes the scenario where multiple cell types have the same eTR/gene pair as significant.[ie classifies as either scenario 4 or 5]

    """
    scenario_same_tr = {}

    # Get the sign of the main cell type's effect size
    main_effect_size_sign = int(np.sign(row_coeff))

    # Loop through all other cell types
    for cell_type in cell_types_with_same_eTR['cell_type'].unique():
        if cell_type == row_cell_type:
            continue  # Skip the main one

        # Get the sign of the secondary cell type's effect size
        secondary_effect_size_sign = int(
            np.sign(
                cell_types_with_same_eTR[cell_types_with_same_eTR['cell_type'] == cell_type]['coeff_meta_fixed'].values[
                    0
                ]
            )
        )

        # Assign scenario based on whether sign matches
        if secondary_effect_size_sign == main_effect_size_sign:
            scenario_same_tr[cell_type] = 4
        else:
            scenario_same_tr[cell_type] = 5

    return scenario_same_tr


def process_same_egene(cell_types_with_same_egene, row_coeff, row_variantid, row_cell_type, ld_path):
    """
    Processes the scenario where multiple cell types have the same gene as signficant (but not necessarily the same TR as being associated with the gene). [ie classifies as either scenario 3/4/5]
    """
    scenario_same_egene = {}

    # Get the sign of the main cell type's effect size
    main_effect_size_sign = int(np.sign(row_coeff))

    # Load the LD matrix for the current gene
    egene_ld = pd.read_csv(f'{ld_path}/{cell_types_with_same_egene["gene_name"].iloc[0]}.csv')

    for cell_type in cell_types_with_same_egene['cell_type'].unique():
        if cell_type == row_cell_type:
            continue  # Skip the main cell type

        # Get the subset for this secondary cell type
        subset = cell_types_with_same_egene[cell_types_with_same_egene['cell_type'] == cell_type]

        # Filter LD matrix to only relevant variantids
        filtered_ld = egene_ld[egene_ld['Unnamed: 0'].isin(subset['variantid'])]

        # Safety check: skip if no matching LD data
        if filtered_ld.empty or row_variantid not in filtered_ld.columns:
            continue

        # Find variant in subset with max LD to row_variantid
        max_ld_value = filtered_ld[row_variantid].max()
        max_row = filtered_ld.loc[filtered_ld[row_variantid].idxmax()]
        max_variantid = max_row['Unnamed: 0']

        if max_ld_value <= 0.5:
            scenario_same_egene[cell_type] = 3  # Not in high LD
        else:
            # Get the sign of the secondary variant's effect
            secondary_coeff = subset[subset['variantid'] == max_variantid]['coeff_meta_fixed'].values
            if secondary_coeff.size == 0:
                continue  # Safety check
            secondary_sign = int(np.sign(secondary_coeff[0]))

            if secondary_sign != main_effect_size_sign:
                scenario_same_egene[cell_type] = 5  # Opposite sign
            else:
                scenario_same_egene[cell_type] = 4  # Same sign

    return scenario_same_egene


def reprocess_scenario_2(
    row_cell_type, scenario_2_celltypes, egene_chrom, row_pos, row_end, egene, row_motif, meta_scen2_path
):
    """
    Reprocesses scenario 2 to check if can be reclassified as scenario 4/5 based on meta-analysis p value, based on the criteria from Cuomo et al. 2025.

    """
    reprocess_scen2_dict = {}

    # Load the meta results for the given cell type
    meta_results = pd.read_csv(f'{meta_scen2_path}/{row_cell_type}/meta_results.tsv', sep='\t')

    for celltype in scenario_2_celltypes:
        metasubset = meta_results[
            (meta_results['chrom'] == egene_chrom)
            & (meta_results['pos'] == row_pos)
            & (meta_results['end'] == row_end)
            & (meta_results['motif'] == row_motif)
            & (meta_results['gene_name'] == egene)
            & (meta_results['cell_type2'] == celltype)
        ]

        if metasubset.empty:
            print(f"[{celltype}] No matching entry found in meta results.")
            continue

        metasubset = metasubset.iloc[0]  # assume unique match

        flag_reclassify = False
        reason = ""

        if metasubset['pval_main'] == 0:
            if metasubset['pval_meta_fixed'] == 0:
                flag_reclassify = True
                reason = "both pval_main and pval_meta_fixed are zero"
        elif metasubset['pval_meta_fixed'] / metasubset['pval_main'] < 1e-5:
            flag_reclassify = True
            reason = f"pval_meta_fixed ({metasubset['pval_meta_fixed']}) is much smaller than pval_main ({metasubset['pval_main']})"

        if flag_reclassify:
            same_sign = (metasubset['coeff_main'] >= 0 and metasubset['coeff_2'] >= 0) or (
                metasubset['coeff_main'] < 0 and metasubset['coeff_2'] < 0
            )
            new_category = 4 if same_sign else 5
            print(
                f"[{celltype}] Reclassified from 2 → {new_category} because {reason}; "
                f"coeff_main = {metasubset['coeff_main']}, coeff_2 = {metasubset['coeff_2']}"
            )
            reprocess_scen2_dict[celltype] = new_category

    return reprocess_scen2_dict


def process_cell_type_specificity(estrs, cell_types, ld_path, meta_scen2_path, associatr_path):
    """
    Helper function to process cell type specificity analysis as per Cuomo et al. 2025 procedure.
    Will call on the other helper functions defined above to process the different scenarios.
    """
    import numpy as np
    import pandas as pd
    from cpg_utils.hail_batch import output_path

    results = []

    for index, row in estrs.iterrows():
        row_dict = {}
        egene = row['gene_name']
        egene_chrom = row['chr']
        row_cell_type = row['cell_type']
        row_pos = row['pos']
        row_end = row['end']
        row_motif = row['motif']
        row_coeff = row['coeff_meta_fixed']
        row_variantid = row['variantid']

        # --- Scenario 4/5: same eTR ---
        cell_types_with_same_eTR = estrs[
            (estrs['gene_name'] == egene)
            & (estrs['pos'] == row_pos)
            & (estrs['end'] == row_end)
            & (estrs['motif'] == row_motif)
        ]
        if len(cell_types_with_same_eTR) > 1:
            same_etr_dict = process_same_tr(cell_types_with_same_eTR, row_cell_type, row_coeff)
            row_dict.update(same_etr_dict)

        # --- Scenario 3/4/5: same gene, different TR ---
        cell_types_with_same_egene = estrs[
            (estrs['gene_name'] == egene) & (~estrs['cell_type'].isin([k for k, v in row_dict.items() if v == 5]))
        ]
        if len(cell_types_with_same_egene) > 1:
            same_egene_dict = process_same_egene(
                cell_types_with_same_egene, row_coeff, row_variantid, row_cell_type, ld_path
            )
            for k, v in same_egene_dict.items():
                if v > row_dict.get(k, -1):
                    row_dict[k] = v

        # --- Scenario 1/2: gene not found in other cell types ---
        cell_types_with_egene = estrs[estrs['gene_name'] == egene]['cell_type']
        cell_types_without_egene = set(cell_types) - set(cell_types_with_egene)
        no_egene_dict = process_scenario_1_2(egene, egene_chrom, cell_types_without_egene, associatr_path)
        for k, v in no_egene_dict.items():
            if v > row_dict.get(k, -1):
                row_dict[k] = v

        # --- Reprocess scenario 2s to 4/5 based on meta pval ---
        keys_with_value_2 = [key for key, value in row_dict.items() if value == 2]
        scen2_dict = reprocess_scenario_2(
            row_cell_type, keys_with_value_2, egene_chrom, row_pos, row_end, egene, row_motif, meta_scen2_path
        )
        for k, v in scen2_dict.items():
            if v > row_dict.get(k, -1):
                row_dict[k] = v

        # --- Final step: make sure all cell types are represented ---
        complete_row = {ct: row_dict.get(ct, np.nan) for ct in cell_types}
        results.append(complete_row)

    # --- Convert final list of dicts into DataFrame ---
    scenario_df = pd.DataFrame(results)
    print(scenario_df)
    scenario_df.columns = cell_types
    scenario_df.index = estrs.index  # to align with the original estrs rows

    # --- Append back to estrs ---
    estrs = pd.concat([estrs, scenario_df.add_prefix("scenario_")], axis=1)

    estrs.to_csv(
        output_path(f'{row_cell_type}/{egene_chrom}/cell_type_specificity_analysis.csv', 'analysis'), index=False
    )


@click.option(
    '--estrs-path',
    type=str,
    required=True,
    help='Path to the estrs file.',
    default='gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv',
)
@click.option(
    '--meta-scen2-path',
    type=str,
    required=True,
    help='Path to the meta scenario 2 file.',
    default='gs://cpg-tenk10k-test-analysis/str/cell-spec-cuomo/meta_results',
)
@click.option(
    '--ld-path',
    type=str,
    required=True,
    help='Path to the LD files directory.',
    default='gs://cpg-tenk10k-test/str/cell-spec-cuomo/v2/ld_files',
)
@click.option(
    '--associatr-path',
    type=str,
    required=True,
    help='Path to the associatr results directory.',
    default='gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2',
)
@click.command()
def main(estrs_path, meta_scen2_path, ld_path, associatr_path):
    estrs = pd.read_csv(estrs_path)
    cell_types = 'CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,CD4_TCM_permuted,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC'
    cell_types = cell_types.split(',')
    estrs['variantid'] = estrs['pos'].astype(str) + estrs['motif']

    b = get_batch(name='Cell-spec-analysis-Cuomo-2025')
    for chrom in range(1, 23):
        for cell_type in cell_types:
            estrs_celltype_chrom = estrs[(estrs['cell_type'] == cell_type) & (estrs['chr'] == f'chr{chrom}')]
            if estrs_celltype_chrom.empty:
                print(f"No data for {cell_type} on chromosome {chrom}, skipping.")
                continue
            if to_path(output_path(f'{cell_type}/{chrom}/cell_type_specificity_analysis.csv', 'analysis')).exists():
                print(f"Cell type specificity analysis for {cell_type} on chromosome {chrom} already exists, skipping.")
                continue
            j = b.new_python_job(
                name=f'Cell spec for {cell_type} and chromosome {chrom}',
            )
            j.call(
                process_cell_type_specificity,
                estrs_celltype_chrom,
                cell_types,
                ld_path,
                meta_scen2_path,
                associatr_path,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()
