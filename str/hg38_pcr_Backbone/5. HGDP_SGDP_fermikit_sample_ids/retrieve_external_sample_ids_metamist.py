import json
import pandas as pd 

from sample_metadata.apis import ParticipantApi,SampleApi

participants = ParticipantApi().get_participants('hgdp')

participant_ids = []
for i in participants:
    participant_ids.append(i['id'])


samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp'], "participant_ids": participant_ids})

external_sample_ids=[] # list of external_sample_ids that are found in cpg hgdp dataset 
for i in samples: 
    external_sample_ids.append(i["external_id"])

sgdp_metadata = open("sgdp_metadata_mappings.csv") #mapping of 279 publicly accessible SGDP IDs to Illumina IDs and Sample IDs

sgdp_metadata = sgdp_metadata.readlines()

sgdp_id_to_external_sample_id={} #dictionary that stores SGDP_ids that have mappings to CPG-HGDP external_sample_ids
sgdp_id_no_mapping = [] #array storing SGDP_ids that have no mapping to metamist cpg-hgdp samples 

for i in sgdp_metadata: 
    if "SGDP_ID" in i: 
        continue # first line, no info 
    sgdp_id, illumina_id, sample_id = i.split(",")
    sgdp_id = sgdp_id.strip()
    illumina_id = illumina_id.strip()
    sample_id = sample_id.strip()
    key = sgdp_id
    if illumina_id in external_sample_ids: #check if SGDP sample is already stored in metamist by its illumina_id attribute
        value = illumina_id
    elif sample_id in external_sample_ids:  #check if SGDP sample is already stored in metamist by its sample_id attribute
        value = sample_id
    else: 
        sgdp_id_no_mapping.append(sgdp_id) #if not stored in metamist, we will append to the array containing SGDP ids that have no mapping to metamist
        continue
    sgdp_id_to_external_sample_id[key] = value 

sgdp_id_no_mapping.append("S_Japanese-1") #SGDP id of LP6005441-DNA_C06 which is listed in metamist but has no crams stored

print(len(sgdp_id_no_mapping)) #148

#print(len(sgdp_id_to_external_sample_id)) #132

#148+132 =280, 1 more than 279 due to the manual addition of LP6005441-DNA_C06 to the sgdp_id_no_mapping array

sgdp_metadata_dataframe = pd.read_csv("sgdp_metadata_mappings.csv")
#subset for SGDP-ids that have no match in metamist
sgdp_metadata_dataframe = sgdp_metadata_dataframe.loc[sgdp_metadata_dataframe['SGDP_ID'].isin(sgdp_id_no_mapping)]
illumina_ids_to_ingest = sgdp_metadata_dataframe[["Illumina_ID"]]

PRJEB9586 = open("filereport_read_run_PRJEB9586_json.txt") #file report of the PRJEB9586 project that contains the 279 SGDP genomes
PRJEB9586 = json.load(PRJEB9586)
sample_accessions_to_ingest=[] #stores sample accession ids that need to be ingested
external_ids_to_ingest=[] #stores sample external ids that correspond to the sample accession ids that need to be ingested 
for i in PRJEB9586: 
    try: 
        link_1,link_2 = i['submitted_ftp'].split(";") 
    except:
        continue #bam file extension
    link_1_sample_id_extension = link_1.split("/")[-1] #eg LP6005441-DNA_C06_1.fastq.gz
    link_1_extension_length = len(link_1_sample_id_extension)
    link_1_sample_id = link_1_sample_id_extension[:link_1_extension_length-11] #removes the '_1.fastq.gz' suffix to obtain the sample id
    if link_1_sample_id in illumina_ids_to_ingest.values: 
        external_ids_to_ingest.append(link_1_sample_id)
        sample_accessions_to_ingest.append(i['sample_accession'])


sample_accessions_to_ingest =list(dict.fromkeys(sample_accessions_to_ingest)) #remove duplicates
external_ids_to_ingest = list(dict.fromkeys(external_ids_to_ingest)) #remove duplicates


print(sample_accessions_to_ingest)
print(len(sample_accessions_to_ingest)) #148
print(external_ids_to_ingest)
print(len(external_ids_to_ingest)) #148


#with open('sgdp_id_to_external_sample_id.txt', 'w') as convert_file:
   #  convert_file.write(json.dumps(sgdp_id_to_external_sample_id))






