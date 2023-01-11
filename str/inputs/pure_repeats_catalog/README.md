# Workflow: 

1) Remove complex repeats from Illumina’s 174k loci catalog (.JSON). Convert .JSON to BED file format. [Script](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/Illumina%20catalog%20to%20BED%20file%20conversion.ipynb) [BED File](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/bed_catalog_without_complex_repeats.bed)

2) Input BED catalog file into BEDTOOLS to extract FASTA sequence of each locus. [Script](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/bed_to_sequence.py) [FASTA](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/catalog_fasta_sequences.fasta.txt) [Hail batch](https://batch.hail.populationgenomics.org.au/batches/420088)

3) Run TandemRepeatFinder v4.09 on the BED catalog (completed locally).

- Motivation: [Broad’s approach](https://gnomad.broadinstitute.org/news/2022-01-the-addition-of-short-tandem-repeat-calls-to-gnomad/), citing improved EH accuracy when looking at only pure repeats  

- Command: `trf catalog_fasta_sequences.fasta.txt 2 7 7 80 10 2 500 -d -h`

- Used default settings apart from `minscore` = 2, to ensure no locus is dropped
[TRF Output file](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/trf_output.dat)

4) Filter TRF output for pure repeats only (100% purity score). [Script](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/Pure%20repeats%20catalog%20.ipynb) [Loci of pure repeats](https://github.com/populationgenomics/sv-workflows/blob/pure_repeats_PR/str/inputs/pure_repeats_catalog/pure_repeats_loci.txt)
