# Workflow

1) Remove complex repeats from Illumina’s 174k loci catalog (.JSON). Convert .JSON to BED file format. 
- Motivation: [Broad’s approach](https://gnomad.broadinstitute.org/news/2022-01-the-addition-of-short-tandem-repeat-calls-to-gnomad/), citing improved EH accuracy when looking at only pure repeats. 
- Output: `bed_catalog_without_complex_repeats.bed`


2) Input BED catalog file into BEDTOOLS to extract FASTA sequence of each locus. - `bed_to_fasta.py`
- Output: `intermediate_files/Illumina_catalog_sequences.fasta.txt`

3) Run `find_pure_repeats.py`; filter for pure loci, filter for loci with motifs between 2-6bp. (164,847 loci are pure)
- Output: multiple files in `intermediate_files`

4) Run `pure_repeats_catalog.py `which identifies if there are additional motifs in the L/R flanks, revises coordinates if necessary, and outputs a BED file with the coordinates of the pure repeats. 
- Output: `intermediate_files/pure_repeat_loci.bed`

5) Run `eh_catalog_creation.py` to create the EH catalog (JSON) based on `intermediate_files/pure_repeat_loci.bed`
- Output: `catalogs/eh_pure_repeats_catalog.json`

6) Run `gangstr_catalog_creation.py` to create the GangSTR catalog (JSON) based on `catalogs/eh_pure_repeats_catalog.json`
- Output: `catalogs/gangstr_pure_repeats_catalog.bed`

