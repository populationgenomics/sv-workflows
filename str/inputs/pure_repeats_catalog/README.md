# Workflow

1) Remove complex repeats from Illumina’s 174k loci catalog (.JSON). Convert .JSON to BED file format. [Script](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/Illumina%20catalog%20to%20BED%20file%20conversion.ipynb) [BED](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/bed_catalog_without_complex_repeats.bed.zip)
  * **Motivation**: [Broad’s approach](https://gnomad.broadinstitute.org/news/2022-01-the-addition-of-short-tandem-repeat-calls-to-gnomad/), citing improved EH accuracy when looking at only pure repeats

2) Input BED catalog file into BEDTOOLS to extract FASTA sequence of each locus. [Script](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/bed_to_fasta.py) [FASTA](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/catalog_fasta_sequences.fasta.txt.zip) [Hail batch](https://batch.hail.populationgenomics.org.au/batches/420088)

3) Run TandemRepeatFinder v4.09 on the FASTA sequences (completed locally).
 * **Command**: `trf catalog_fasta_sequences.fasta.txt 2 7 7 80 10 2 500 -d -h`
 * (Used default settings apart from `minscore` = 2, to ensure no locus is dropped
[TRF Output file](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/trf_output.dat.zip))

4) Filter TRF output for pure repeats only (100% purity score). [Script](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/pure_repeats_script.py) [Loci of pure repeats](https://github.com/populationgenomics/sv-workflows/blob/main/str/inputs/pure_repeats_catalog/pure_repeat_loci.txt.zip)
