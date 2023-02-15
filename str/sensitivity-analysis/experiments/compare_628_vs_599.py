pemberton = open("Pemberton_AdditionalFile1_11242009.txt")
hipstr_regions = open("marshfield_regions.bed")
pemberton = pemberton.readlines()
hipstr_regions = hipstr_regions.readlines()

pemberton_loci=[]
for i in pemberton: 
    pemberton_loci.append(i.split()[3])

print(len(pemberton_loci)) #628, representing 627 loci and one header 

hipstr_loci = []
for j in hipstr_regions: 
    hipstr_loci.append(j.split()[5])
print(len(hipstr_loci))
#print(hipstr_loci[5])

excluded_loci = []
for k in pemberton_loci: 
    if k not in hipstr_loci: 
        excluded_loci.append(k)

print((excluded_loci))