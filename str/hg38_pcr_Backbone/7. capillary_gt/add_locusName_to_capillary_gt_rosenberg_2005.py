import pandas as pd 

capillary = pd.read_csv("capillary_genotypes_rosenberg2005.csv")
print(capillary)

alternate_name = pd.read_csv("Pemberton_AdditionalFile1_11242009.txt", sep = "\t")
alternate_name = alternate_name[["locusName", "nameInGenotypeDataFile", "alternateName"]]
capillary = pd.merge(capillary, alternate_name, left_on = "locus_id", right_on = "nameInGenotypeDataFile")
print(capillary)




capillary.to_csv("capillary_genotypes_with_locusNames_rosenberg2005.csv")