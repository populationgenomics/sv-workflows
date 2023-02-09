"""
This script filters for HGDP ids that fit the French, Han or Yoruba population groups. 
"""

from sample_metadata.apis import AnalysisApi, ParticipantApi

data = ParticipantApi().get_participants('hgdp')
population_groups = ["French", "Yoruba", "Han", "Northern Han"]
sample_ids = []


for i in data: 
    if i["meta"]["Population name"] in population_groups: 
        sample_ids.append(i['external_id'])
        print(i["meta"]["Population name"])

print(len(sample_ids))
print(sample_ids)
