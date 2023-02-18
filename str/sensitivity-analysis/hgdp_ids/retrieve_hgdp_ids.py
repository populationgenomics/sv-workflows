"""
This script filters for HGDP ids that fit the French, Han or Yoruba population groups. 
"""

import sample_metadata
from sample_metadata.apis import SeqrApi, ParticipantApi,SampleApi
external_participant_ids=["HGDP00019","HGDP00027","HGDP00058","HGDP00090","HGDP00124","HGDP00125","HGDP00157","HGDP00160","HGDP00195","HGDP00208","HGDP00216","HGDP00232","HGDP00286","HGDP00328","HGDP00338","HGDP00428","HGDP00449","HGDP00457","HGDP00461","HGDP00474","HGDP00476","HGDP00526","HGDP00530","HGDP00540","HGDP00541","HGDP00543","HGDP00545","HGDP00547","HGDP00548","HGDP00549","HGDP00550","HGDP00551","HGDP00552","HGDP00553","HGDP00554","HGDP00555","HGDP00556","HGDP00569","HGDP00597","HGDP00616","HGDP00650","HGDP00656","HGDP00660","HGDP00702","HGDP00706","HGDP00713","HGDP00717","HGDP00722","HGDP00725","HGDP00737","HGDP00749","HGDP00773","HGDP00783","HGDP00785","HGDP00796","HGDP00798","HGDP00846","HGDP00852","HGDP00855","HGDP00857","HGDP00887","HGDP00903","HGDP00915","HGDP00928","HGDP00932","HGDP00951","HGDP00956","HGDP00987","HGDP00991","HGDP01012","HGDP01018","HGDP01028","HGDP01030","HGDP01032","HGDP01034","HGDP01035","HGDP01044","HGDP01047","HGDP01078","HGDP01079","HGDP01095","HGDP01098","HGDP01153","HGDP01163","HGDP01168","HGDP01172","HGDP01179","HGDP01188","HGDP01191","HGDP01198","HGDP01199","HGDP01203","HGDP01211","HGDP01215","HGDP01223","HGDP01228","HGDP01240","HGDP01242","HGDP01246","HGDP01250","HGDP01253","HGDP01274","HGDP01297","HGDP01306","HGDP01312","HGDP01314","HGDP01315","HGDP01320","HGDP01323","HGDP01333","HGDP01335","HGDP01338","HGDP01345","HGDP01350","HGDP01355","HGDP01364","HGDP01365","HGDP01401","HGDP01402","HGDP01414","HGDP01417"]

data = ParticipantApi().get_participants('hgdp', body_get_participants={'external_participant_ids':external_participant_ids})
#population_groups = ["French", "Yoruba", "Han", "Northern Han"]

participants_id =[]
#for i in data: 
 #   if i["meta"]["Population name"] in population_groups: 
     #   participant_ids.append(i['id']) 
for i in data:
    participants_id.append(i['id'])


samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp'], "participant_ids": participants_id})
print(samples[0])

sample_ids=[]
for i in samples: 
    sample_ids.append(i["external_id"])



print(len(sample_ids))
print(sample_ids)

