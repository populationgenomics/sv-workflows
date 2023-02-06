
import re

marshfield_regions = open("marshfield_regions.bed")
marshfield_regions = marshfield_regions.readlines()
marshfield_dict={}
for j in marshfield_regions: #599 loci
    attributes = j.split()
    locus = attributes[5]
    chr = "chr"+attributes[0]
    start_coord = attributes[1]
    end_coord = attributes[2]
    marshfield_dict[locus] = (chr,start_coord, end_coord)

print(len(marshfield_dict)) #599

file = open("Pemberton_AdditionalFile1_11242009.txt")
hg19_coordinates = open("hg19_coordinates.bed", "w")

file = file.readlines()
print(len(file)) #628 

pure_repeat_locus_motif_dict ={} 
counter=0
for i in file: 
    attributes = i.split()
    locus = attributes[3]
    motif = attributes[9]
    if motif.count("(") ==1: #indicative of a pure simple repeat 
        motif = re.sub(r'[^a-zA-Z]', '', motif)
        pure_repeat_locus_motif_dict[locus] = motif 

print(len(pure_repeat_locus_motif_dict)) #num pure simple repeats (488)


for locus in pure_repeat_locus_motif_dict: 
    if locus in marshfield_dict.keys():
        chr, start, end = marshfield_dict[locus]
        hg19_coordinates.write(chr+"\t"+start+"\t"+end+"\t"+locus+"\n")
        counter+=1
    else:
        if locus == "SRA": #must be done manually because not found in marshfield_dict
            hg19_coordinates.write("chr2\t1493422\t1493453\tSRA\n")
            counter+=1
        print(locus) # loci without matching coordinates in Marshfield panel 
print(counter) #482 loci (ie not all the loci are found in the Marshfield coordinates; 1 locus was manually placed; other 6 locus failed to be uniquely mapped to ref genome using UCSC in silico CPR


hg19_coordinates.close()

hg38_coordinates= open("hg38_coordinates.bed") #converted hg19_coordinates.bed to hg38 (total 482 loci)
hg38_coordinates = hg38_coordinates.readlines() 

hg38_dict = {}
for k in hg38_coordinates:
    chr, start, end, locus, num = k.split() 
    hg38_dict[locus] = (chr, start, end)

print(len(hg38_dict)) #482


trf_output = open("hg38_coordinates.fasta.2.5.7.80.10.16.2000.dat")
trf_output = trf_output.readlines()
n = 0 
trf_dict={}
while n< len(trf_output):
    if "Sequence" in trf_output[n]:
        locus = trf_output[n][10:].rstrip()
        n+=7 
        trf_motif = trf_output[n].split()[13]
        trf_dict[locus] = trf_motif
    n+=1
print(trf_dict["chr22:47865081-47865104"])
print(len(trf_dict))#482 loci

gangstr_catalog = open("marshfield_regions_gangstr_untrimmed_coordinates.bed", "w")
counter_2=0
for locus in hg38_dict:
    chr, start, end = hg38_dict[locus]
    trf_dict_key = f'{chr}:{start}-{end}'
    start = str(int(start)+1) #Gangstr is 1-based
    motif = trf_dict[trf_dict_key]
    motif_len = str(len(motif))
    gangstr_catalog.write(chr+"\t"+start+"\t"+end+"\t"+motif_len+"\t"+motif+"\n")
    counter_2+=1
print(counter_2) #482 loci
gangstr_catalog.close()


