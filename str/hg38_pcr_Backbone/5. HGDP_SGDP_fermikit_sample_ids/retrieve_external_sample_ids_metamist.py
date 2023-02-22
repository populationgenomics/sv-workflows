import json

from sample_metadata.apis import SeqrApi, ParticipantApi,SampleApi

participants = ParticipantApi().get_participants('hgdp')

participant_ids = []
for i in participants:
    participant_ids.append(i['id'])


samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp'], "participant_ids": participant_ids})
print(samples[0])

external_sample_ids=[]
for i in samples: 
    external_sample_ids.append(i["external_id"])

fermikit_mappings = open("Fermikit_ID_mapping.csv") 

fermikit_mappings = fermikit_mappings.readlines()

sgdp_id_to_external_sample_id={}

for i in fermikit_mappings: 
    if "SGDP_ID" in i: 
        continue # first line, no info 
    sgdp_id, illumina_id, sample_id = i.split(",")
    sgdp_id = sgdp_id.strip()
    illumina_id = illumina_id.strip()
    sample_id = sample_id.strip()
    key = sgdp_id
    if illumina_id in external_sample_ids: 
        value = illumina_id
    elif sample_id in external_sample_ids: 
        value = sample_id
    else: 
        value = "no mapping" 
        continue
    sgdp_id_to_external_sample_id[key] = value 

print(len(sgdp_id_to_external_sample_id))#121
print(sgdp_id_to_external_sample_id.values()) #sample IDs to run on GangSTR and EH 

with open('sgdp_id_to_external_sample_id.txt', 'w') as convert_file:
     convert_file.write(json.dumps(sgdp_id_to_external_sample_id))






