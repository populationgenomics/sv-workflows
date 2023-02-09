"""
This script filters for HGDP ids that fit the French, Han or Yoruba population groups. 
"""
from sample_metadata.apis import SampleApi, ParticipantApi

data = ParticipantApi().get_participants('hgdp-test')
print(data[0])
population_groups = ["French", "Yoruba", "Han", "Northern Han"]
participant_ids = []


for i in data: 
    if i["meta"]["Population name"] in population_groups: 
        participant_ids.append(i['id']) 
sample_ids=[]
samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp-test'], "participant_ids": participant_ids})
print(samples[1])
for i in samples: 
    sample_ids.append(i["id"])



print(len(sample_ids))
print(sample_ids)
