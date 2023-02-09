"""
This script filters for HGDP ids that fit the French, Han or Yoruba population groups. 
"""
import sample_metadata
from sample_metadata.apis import SampleApi, ParticipantApi

data = ParticipantApi().get_participants('hgdp')
print(data[0])
population_groups = ["French", "Yoruba", "Han", "Northern Han"]
participant_ids = []


for i in data: 
    if i["meta"]["Population name"] in population_groups: 
        participant_ids.append(i['id'])

body_get_samples = SampleApi().BodyGetSamples(
        sample_ids=[
            "",
        ],
        meta={},
        participant_ids= participant_ids,
        project_ids=[
            "hgdp",
        ],
        active=True,
    ) 

SampleApi().get_samples(body_get_samples={'project_ids':['hgdp'], "participant_ids": participant_ids})



print(len(sample_ids))
print(sample_ids)
