"""
This script identifies if there are additional repeats in the flanks of pure repeat loci, and modifies the coordinates to include them as necessary. 

Required packages: str_analysis
python3 -m pip install --upgrade str_analysis
"""

from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence

## File loading of catalog with L/R flanks that are one motif length long. 
d= open('intermediate_files/catalog_with_flanks_one_motif_length.fasta')
flank_catalog_one_motif_length = d.readlines()
print(len(flank_catalog_one_motif_length))#989094 lines (ie each of the 164849 loci gets 6 lines)

## Dictionary creation
flank_catalog_one_motif_length_dict={}
i=0
while i<len(flank_catalog_one_motif_length)-5:
    if ">" in flank_catalog_one_motif_length[i]:
        key = flank_catalog_one_motif_length[i][1:].rstrip() # removes the ">"
        value = (flank_catalog_one_motif_length[i+3].rstrip(),flank_catalog_one_motif_length[i+5].rstrip()) # tuple (FASTA + R flank, FASTA + L flank)
        flank_catalog_one_motif_length_dict[key] = value 
    i+=6
print(len(flank_catalog_one_motif_length_dict)) #164849

## Pure repeat catalog dictionary for comparison
f = open("intermediate_files/pure_repeat_catalog_not_final.txt")
pure_catalog = f.readlines()
j=0
pure_catalog_dict={}
while j<len(pure_catalog):
    entry = pure_catalog[j].rstrip()
    locus, motif, repeat_count = entry.split(" ")
    pure_catalog_dict[locus] = (motif, repeat_count)
    j+=1

## Look into L/R flanks for additional repeats 
left_flank_one_motif_added_repeat_loci=[]
right_flank_one_motif_added_repeat_loci =[]
both_flanks_one_motif_added_repeat_loci=[]

for locus in flank_catalog_one_motif_length_dict:
    right_flank, left_flank = flank_catalog_one_motif_length_dict[locus]
    motif, num_repeats = pure_catalog_dict[locus]
    right_flank_num_repeats= extend_repeat_into_sequence(motif, right_flank)[0] #returns the num_repeats of the longest stretch of repeats with the motif. Sequence input must start with the motif, else returns 0. 
    left_flank_num_repeats = extend_repeat_into_sequence(motif, left_flank)[0]
    if (left_flank_num_repeats> int(num_repeats) and right_flank_num_repeats> int(num_repeats) ):
        both_flanks_one_motif_added_repeat_loci.append(locus)
    elif (left_flank_num_repeats> int(num_repeats)):
        left_flank_one_motif_added_repeat_loci.append(locus)
    elif (right_flank_num_repeats> int(num_repeats)):
        right_flank_one_motif_added_repeat_loci.append(locus)

print(len(left_flank_one_motif_added_repeat_loci)) #31
print(len(right_flank_one_motif_added_repeat_loci)) #0
print(len(both_flanks_one_motif_added_repeat_loci)) #0

# optional 
left_flank_added_repeat_loci_file = open("intermediate_files/left_flank_one_motif_added_repeat_loci.txt", "w")
left_flank_added_repeat_loci_file.write(str(left_flank_one_motif_added_repeat_loci))
left_flank_added_repeat_loci_file.close()

#right_flank_added_repeat_loci and both_flanks_added_repeat_loci are empty 

## Because there were some additional repeats detected in the left flanks, will consider L/R flanks that are 2 motifs long now. 

## File loading of catalog with L/R flanks that are two motif lengths long. 
g= open('intermediate_files/catalog_with_flanks_two_motif_lengths.fasta')
flank_catalog_two_motif_lengths = g.readlines()
print(len(flank_catalog_two_motif_lengths))#989094 lines (ie each of the 164849 loci gets 6 lines)

## Dictionary creation of catalog with L/R flanks that are two motif lengths long: 
flank_catalog_two_motif_lengths_dict={} #164849

i=0
while i<len(flank_catalog_two_motif_lengths)-5:
    if ">" in flank_catalog_two_motif_lengths[i]:
        key = flank_catalog_two_motif_lengths[i][1:].rstrip() #removes the ">"
        value = (flank_catalog_two_motif_lengths[i+3].rstrip(),flank_catalog_two_motif_lengths[i+5].rstrip()) # tuple (FASTA + R flank, FASTA + L flank)
        flank_catalog_two_motif_lengths_dict[key] = value 
    i+=6
print(len(flank_catalog_two_motif_lengths_dict))

## Look into L flank for additional repeats 
left_flank_two_motif_added_repeat_loci=[]

for locus in flank_catalog_two_motif_lengths_dict:
    right_flank, left_flank = flank_catalog_two_motif_lengths_dict[locus]
    motif, num_repeats = pure_catalog_dict[locus]
    left_flank_num_repeats = extend_repeat_into_sequence(motif, left_flank)[0]
    if (left_flank_num_repeats> int(num_repeats)): #just checking L flank now because the R flank yielded no additional repeats in prev. iteration
        left_flank_two_motif_added_repeat_loci.append(locus)
print(len(left_flank_two_motif_added_repeat_loci)) #0, suggesting that all the repeats in the left flanks are only one copy

## Edit pure repeats catalog to append repeats in L flank (31 loci)
revised_pure_catalog_dict ={}

for locus in pure_catalog_dict:
    motif, repeat_count = pure_catalog_dict[locus]
    if locus in left_flank_one_motif_added_repeat_loci:
        motif_length = len(motif)
        chr = locus.split(":")[0]
        coordinates = locus.split(":")[1]
        start_coordinate = int(coordinates.split("-")[0])
        start_coordinate_minus_motif_length = int(coordinates.split("-")[0]) -motif_length
        end_coordinate = int(coordinates.split("-")[1])
        revised_locus = chr+":"+str(start_coordinate_minus_motif_length)+"-"+str(end_coordinate)
        revised_repeat_count = int(repeat_count)+1
        revised_pure_catalog_dict[revised_locus] = (motif,revised_repeat_count)
    else:
        revised_pure_catalog_dict[locus]= (motif, repeat_count)
print(len(revised_pure_catalog_dict)) #164849

#testing out revised pure catalog: 
#print(pure_catalog_dict["chr10:10022442-10022450"]) #('CTGC', '2')
#print(revised_pure_catalog_dict["chr10:10022438-10022450"]) #('CTGC', 3)

#print(pure_catalog_dict["chr10:129364172-129364184"]) #('AATA', '3')
#print(revised_pure_catalog_dict["chr10:129364168-129364184"]) #('AATA', 4)

# Write out revised pure repeat catalog as a BED file in the format chr \t start_coordinate \t end_coordinate \t motif_length \t motif \t num_repeats 
pure_repeats_bed_catalog = open("intermediate_files/pure_repeats_loci.bed", "w")
for p in revised_pure_catalog_dict:
    motif, repeat_count = revised_pure_catalog_dict[p]
    motif_length = len(motif)
    chr = p.split(":")[0]
    coordinates = p.split(":")[1]
    start_coordinate = coordinates.split("-")[0]
    end_coordinate = coordinates.split("-")[1]
    pure_repeats_bed_catalog.write(chr+"\t" + str(start_coordinate) +"\t"+str(end_coordinate)+"\t"+str(motif_length)+"\t"+str(motif)+"\t"+str(repeat_count)+"\n") 
pure_repeats_bed_catalog.close()


