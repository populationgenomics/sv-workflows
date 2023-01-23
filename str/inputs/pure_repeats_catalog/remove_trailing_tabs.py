f = open("catalogs/gangstr_pure_repeats_catalog.bed")
file= f.readlines()

print(file[1].split("\t")[0])



with open("gangstr_bed_intersect.bed", "w") as d:

    for i in file: 
        chr = i.split("\t")[0]
        start = i.split("\t")[1]
        start = str(int(start)-1) #make 1-based coordinates back to 0-based for gencode
        end = i.split("\t")[2]
        d.write(chr+"\t"+start+"\t"+end+"\n")
d.close()
