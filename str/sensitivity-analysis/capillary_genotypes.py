import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd 

file = open("combinedmicrosats_627loci_1048indivs_numRpts.stru.txt")
file = file.readlines()

loci = file[0].split()
print(len(loci)) #627

motif_length = file[1].split()
print(len(motif_length)) #627

num_str_regions = file[2].split()
print(len(num_str_regions)) #627 

repeat_regularity = file[3].split()

loci_dict={}
i = 0
while i<len(loci):
    loci_dict[ loci[i]] = (motif_length[i], num_str_regions[i],repeat_regularity[i])
    i+=1

j=4 #start of genotype read outs 
df_main = pd.DataFrame(columns = ['sample_id', 'locus_id', 'genotype', 'population_id', 'population_name', 'geographic_population', 'predefined_region', 'regularity', 'motif_length'])

while j<len(file): 
    attributes = file[j].split()
    sample_id = attributes[0]
    population_id = attributes[1]
    population_name = attributes[2]
    geographic_population = attributes[3]
    predefined_region = attributes[4]
    k = 0 
    while k<627:
        motif_length, num_str_regions,repeat_regularity = loci_dict[loci[k]]
        df_temp = pd.DataFrame({
            'sample_id': [sample_id],
            'locus_id': [loci[k]],
            'genotype':[attributes[5+k]],
            'population_id': [population_id],
            'population_name': [population_name],
            'geographic_population': [geographic_population],
            'predefined_region': [predefined_region] ,
            'regularity': [repeat_regularity],
            'motif_length':[motif_length]
         })
        df_main = df_main.append(df_temp)
        k+=1
    j+=1
    print(j)
df_main.to_csv("capillary_genotypes.csv", index = False)




