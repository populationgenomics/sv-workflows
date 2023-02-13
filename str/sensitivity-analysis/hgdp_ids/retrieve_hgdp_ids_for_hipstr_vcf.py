import sample_metadata
from sample_metadata.apis import SeqrApi, ParticipantApi,SampleApi

file = open("str_sensitivity-analysis_hipstr_untrimmed_coordinates_1_based_91_hgdp_genomes_hipstr.vcf")
writer = open("str_sensitivity-analysis_hipstr_untrimmed_coordinates_1_based_91_hgdp_genomes_hipstr_with_cpgids.csv", "w")
file = file.readline()


attributes = file.split("\t")
first_line =[]
j=0
while j<len(attributes):
    if j<9: 
        first_line.append(attributes[j])
    else: 
        external_id = attributes[j].rstrip()
        cpg_id = SampleApi().get_sample_by_external_id(project ='hgdp-test', external_id = external_id)['id']
        first_line.append(cpg_id)
    j+=1
for k in first_line: 
    writer.write(k+"\t")
writer.write("\n")
        

