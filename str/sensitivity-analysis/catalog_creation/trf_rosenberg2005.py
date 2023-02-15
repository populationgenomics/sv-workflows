trf_output = open("hg38_coordinates.fasta.2.5.7.80.10.16.2000.dat")
trf_output = trf_output.readlines()
n = 0 
trf_dict={}
while n< len(trf_output):
    if "Sequence" in trf_output[n]:
        locus = trf_output[n][10:].rstrip()
        n+=7 
        trf_motif = trf_output[n].split()[13]
        trf_repeat_seq = trf_output[n].split()[14]
        trf_dict[locus] = (trf_motif, trf_repeat_seq)
    n+=1
print(trf_dict["chr22:47865081-47865104"])
print(len(trf_dict))#482 loci

writer = open("trf_output_for_rosenberg2005.csv", "w")
writer.write("locus,motif,repeat_sequence\n")
for i in trf_dict: 
    motif, repeat_seq = trf_dict[i]
    writer.write(f'{i},{motif},{repeat_seq}\n')

