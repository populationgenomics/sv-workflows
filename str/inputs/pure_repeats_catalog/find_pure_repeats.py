"""
This script identifies pure repeat loci in the Illumina catalog. It then creates left and right flanking sequences of pure repeat loci, which will be used as input
in `pure_repeats_catalog.py` to see if there are additional repeats in the flanks.

Required packages: str_analysis
python3 -m pip install --upgrade str_analysis
"""

from str_analysis.utils.find_repeat_unit import find_repeat_unit_kmer

## File loading 
f = open('intermediate_files/Illumina_catalog_sequences.fasta.txt') #original Illumina catalog converted to FASTA sequences of each coordinate pair; removed prior 6 loci containing compound repeats 
orig_catalog = f.readlines()

## Illumina catalog dictionary creation 
orig_catalog_dict ={}
i=0
while i<len(orig_catalog):
    if ">" in orig_catalog[i]:
        key = orig_catalog[i][1:].rstrip()#removes the ">"
        value = orig_catalog[i+1].rstrip()
        orig_catalog_dict[key] = value 
    i+=1
print(len(orig_catalog_dict)) #174287

## Replace value of the dictionary with output from find_repeat_unit_kmer()
for j in orig_catalog_dict:
    kmer_tuple = find_repeat_unit_kmer(orig_catalog_dict[j]) #tuple is in the form (motif, number of repeats found in the sequence)
    orig_catalog_dict[j] = kmer_tuple

## Find loci to exclude from the catalog - loci wtih impure repeat sequences, loci made up of homopolymers
excluded_loci = []
for k in orig_catalog_dict:
    repeat_unit, repeat_count = orig_catalog_dict[k]
    if repeat_count == 1:#ie impure repeat sequence 
        excluded_loci.append(k)
    elif len(repeat_unit)==1:#homopolymers
        excluded_loci.append(k)
    elif len(repeat_unit)>6: #remove loci with motifs greater than 6bp (STR definition is 2-6bp motif)
        excluded_loci.append(k)

## Write out excluded_loci into an output file 
excluded_loci_catalog = open("intermediate_files/excluded_loci_catalog.txt","w")
for m in excluded_loci: 
    motif, repeat_count = orig_catalog_dict[m]
    excluded_loci_catalog.write(m + " "+motif+" "+str(repeat_count)+"\n")
excluded_loci_catalog.close()
print(len(excluded_loci)) #9438

## Remove the excluded_loci from the Illumina catalog dictionary 
for m in excluded_loci:
    del orig_catalog_dict[m]

print(len(orig_catalog_dict)) #164849 pure repeat loci 

## Write out the pure repeat loci to an output file 
pure_repeat_catalog = open("intermediate_files/pure_repeat_catalog_not_final.txt","w")
for n in orig_catalog_dict:
     motif, repeat_count = orig_catalog_dict[n]
     pure_repeat_catalog.write(n +" "+motif + " "+ str(repeat_count)+"\n")
pure_repeat_catalog.close()
   
## Write out pure repeat loci with L/R flanks (size = motif length of STR identified at the locus) in BED format 
bed_catalog_v1 = open("intermediate_files/catalog_with_flanks_one_motif_length.bed", "w")
for o in orig_catalog_dict:
    motif, repeat_count = orig_catalog_dict[o]
    motif_length = len(motif)
    chr = o.split(":")[0]
    coordinates = o.split(":")[1]
    start_coordinate = int(coordinates.split("-")[0])
    start_coordinate_minus_motif_length = int(coordinates.split("-")[0]) -motif_length
    end_coordinate = int(coordinates.split("-")[1])
    end_coordinate_plus_motif_length = end_coordinate + motif_length
    bed_catalog_v1.write(chr+"\t" + str(start_coordinate) +"\t"+str(end_coordinate)+"\n") #original FASTA sequence
    bed_catalog_v1.write(chr+"\t" + str(start_coordinate) +"\t"+str(end_coordinate_plus_motif_length)+"\n") #FASTA sequence + R flank
    bed_catalog_v1.write(chr+"\t" + str(start_coordinate_minus_motif_length) +"\t"+str(end_coordinate)+"\n") # FASTA sequence + L flank
bed_catalog_v1.close()

## Write out pure repeat loci with L/R flanks (size = 2* motif length of STR identified at locus) in BED format 
bed_catalog_v2 = open("intermediate_files/catalog_with_flanks_two_motif_lengths.bed", "w")
for p in orig_catalog_dict:
    motif, repeat_count = orig_catalog_dict[p]
    motif_length = len(motif)
    chr = p.split(":")[0]
    coordinates = p.split(":")[1]
    start_coordinate = int(coordinates.split("-")[0])
    start_coordinate_minus_two_motif_lengths = int(coordinates.split("-")[0]) -2*motif_length
    end_coordinate = int(coordinates.split("-")[1])
    end_coordinate_plus_two_motif_lengths = end_coordinate + 2*motif_length
    bed_catalog_v2.write(chr+"\t" + str(start_coordinate) +"\t"+str(end_coordinate)+"\n") #original FASTA sequence 
    bed_catalog_v2.write(chr+"\t" + str(start_coordinate) +"\t"+str(end_coordinate_plus_two_motif_lengths)+"\n") #FASTA sequence + R flank
    bed_catalog_v2.write(chr+"\t" + str(start_coordinate_minus_two_motif_lengths) +"\t"+str(end_coordinate)+"\n") #FASTA sequence + L flank 
bed_catalog_v2.close()






