import re
import json
regex = re.compile('[^a-zA-Z]')
pemberton = open("Pemberton_AdditionalFile1_11242009.txt")
pemberton = pemberton.readlines() 
pemberton_dict = {}
i = 1 #ignore header line
while i<len(pemberton): 
    attributes = pemberton[i].split("\t")
    alternateName = attributes[3]
    repeatStructure = regex.sub('', attributes[9])
    pemberton_dict[alternateName] = repeatStructure
    i+=1

trf_output = open("trf_input.txt.2.30.30.80.10.16.2000.dat")
trf_output = trf_output.readlines()
n = 0 
trf_dict={}
while n< len(trf_output):
    if "Sequence" in trf_output[n]:
        locus = (trf_output[n][10:].rstrip()).split("_")[1]
        n+=7 
        trf_motif = trf_output[n].split()[13]
        if trf_motif == pemberton_dict[locus]: 
            trf_dict[locus] = trf_output[n].rstrip()
        else: 
            while (n<len(trf_output)-1):
                n+=1
                try:
                    trf_motif = trf_output[n].split()[13]
                    if trf_motif == pemberton_dict[locus]: 
                        trf_dict[locus] = trf_output[n].rstrip()
                        break
                except: 
                     trf_dict[locus] = "no matching motif"
                     break

    n+=1
with open('filtered_trf_output.txt', 'w') as convert_file:
     convert_file.write(json.dumps(trf_dict))