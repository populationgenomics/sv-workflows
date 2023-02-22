
trf_output = open("trf_input_antisense.txt.2.7.7.80.10.8.500.dat")
trf_output = trf_output.readlines()
n = 0 
trf_dict={}
while n< len(trf_output):
    if "Sequence" in trf_output[n]:
        locus = (trf_output[n][10:].rstrip())
        n+=7 
        try:
            trf_motif = trf_output[n].split()[13]
        except: 
            print(trf_output[n-7])
            continue
        trf_dict[locus] = trf_motif
    n+=1
with open('trf_output_antisense.txt', 'w') as convert_file:
    convert_file.write("sequence,motif\n")
    for i in trf_dict: 
        convert_file.write(i + "," + trf_dict[i]+"\n")
    convert_file.write("chr7:131175019-131175027,ATAG\n")

# chr7:131175019-131175027 cannot be added because the sequence is too short, must add manually 