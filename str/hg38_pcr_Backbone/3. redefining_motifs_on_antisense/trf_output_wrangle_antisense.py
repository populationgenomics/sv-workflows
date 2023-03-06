
trf_output = open("trf_input_antisense_debugged.txt.2.7.7.80.10.8.500.dat")
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
with open('trf_output_antisense_debugged.txt', 'w') as convert_file:
    convert_file.write("sequence,motif\n")
    for i in trf_dict: 
        convert_file.write(i + "," + trf_dict[i]+"\n")