import pandas as pd 

capillary = pd.read_csv("capillary_genotypes_rosenberg2005.csv")
print(capillary)

alternate_name = pd.read_csv("Pemberton_AdditionalFile1_11242009.txt", sep = "\t")
alternate_name = alternate_name[["locusName", "nameInGenotypeDataFile", "alternateName"]]
capillary = pd.merge(capillary, alternate_name, left_on = "locus_id", right_on = "nameInGenotypeDataFile")
print(capillary)
"""
coordinates = pd.read_csv("hg38_coordinates.bed", sep = "\t", names = ["chr", "start", "end", "locus_id", "num"])
print(coordinates)
coordinates['start'] = (coordinates['start'].astype(int)).astype(str)
coordinates['end'] = coordinates['end'].astype(int).astype(str)
capillary = pd.merge(coordinates, capillary,left_on = "locus_id", right_on = "alternateName" ,how = "left")
print(capillary)
capillary= capillary.drop(columns=["locus_id_x", "locus_id_y", "num"])
"""




capillary.to_csv("capillary_genotypes_with_locusNames_rosenberg2005.csv")