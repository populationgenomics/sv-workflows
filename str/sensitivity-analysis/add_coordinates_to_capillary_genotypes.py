import pandas as pd 

capillary = pd.read_csv("capillary_genotypes.csv")
#print(capillary)

alternate_name = pd.read_csv("Pemberton_AdditionalFile1_11242009.txt", sep = "\t")
alternate_name = alternate_name[["locusName", "alternateName"]]
capillary = pd.merge(capillary, alternate_name, left_on = "locus_id", right_on = "locusName")
#print(capillary)

coordinates = pd.read_csv("hg38_coordinates.bed", sep = "\t", names = ["chr", "start", "end", "locus_id", "num"])

capillary = pd.merge(capillary, coordinates, left_on = "alternateName", right_on = "locus_id" ,how = "left")

capillary= capillary.drop(columns=["locus_id_x", "locus_id_y", "num"])

capillary.to_csv("capillary_genotypes_with_coordinates.csv")