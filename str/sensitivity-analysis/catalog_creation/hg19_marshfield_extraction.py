marshfield_coordinates = open("hg19_marshfield_coordinates_only.bed", "w")
marshfield_regions = open("marshfield_regions.bed")
marshfield_regions = marshfield_regions.readlines()
counter = 0
for j in marshfield_regions: 
    attributes = j.split()
    locus = attributes[5]
    chr = "chr"+attributes[0]
    start_coord = attributes[1]
    end_coord = attributes[2]
    marshfield_coordinates.write(chr+"\t"+start_coord+"\t"+end_coord+"\n")
    counter+=1
print(counter)