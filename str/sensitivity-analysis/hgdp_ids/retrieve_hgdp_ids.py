"""
This script filters for HGDP ids that fit the French, Han or Yoruba population groups. 
"""
import sample_metadata
from sample_metadata.apis import SeqrApi, ParticipantApi
"""
data = ParticipantApi().get_participants('hgdp')
print(data[0])
population_groups = ["French", "Yoruba", "Han", "Northern Han"]
participant_ids = []


for i in data: 
    if i["meta"]["Population name"] in population_groups: 
        participant_ids.append(i['id']) 
sample_ids=[]
samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp'], "participant_ids": participant_ids})
for i in samples: 
    sample_ids.append(i["external_id"])



print(len(sample_ids))
print(sample_ids)
"""

datum = SeqrApi().get_external_participant_id_to_internal_sample_id('hgdp-test')
datum = dict(datum)
print(datum["HGDP00511"])
