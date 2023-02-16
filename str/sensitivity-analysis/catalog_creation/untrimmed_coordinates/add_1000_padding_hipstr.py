file = open("catalogs/hg38_hipstr_catalog_untrimmed_coordinates_1_based.bed")
file = file.readlines()

padded_catalog = open("catalogs/hg38_hipstr_catalog_untrimmed_coordinates_1_based_3_col.bed", "w")

for i in file: 
    attributes = i.split(" ") 
    start_coord = int(attributes[1])
    assert start_coord >0
    end_coord = int(attributes[2])
    padded_catalog.write(attributes[0]+"\t" +str(start_coord)+"\t"+ str(end_coord)+"\n")


