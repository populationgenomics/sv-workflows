import re
import json
regex = re.compile('[^a-zA-Z]')
pemberton = open("Pemberton_AdditionalFile1_11242009.txt")
pemberton = pemberton.readlines() 
pemberton_dict = {}

def get_metrics_match(metric_list: list[str]) -> str:
    """
    takes list of metrics line and gives the line with the matching motif to Pemberton_dict AND longest repeat sequence if multiple metric lines
    have the same matching motif
    """
    metrics_with_matching_motif= []
    for each_metric in metric_list:
        trf_motif = each_metric.split()[13]
        if trf_motif == pemberton_dict[locus]: 
            metrics_with_matching_motif.append(each_metric.rstrip())
    if len(metrics_with_matching_motif) ==0: 
        return "no matching motif"
    if len(metrics_with_matching_motif) ==1: 
        return metrics_with_matching_motif[0]
    metric_data=[]
    for each_metric in metrics_with_matching_motif: 
        split_line = each_metric.rstrip().split()
        metric_data.append((each_metric, float(split_line[3])))
    max_metrics, num_repeats = max(metric_data, key =lambda x:x[1])
    return max_metrics


i = 1 #ignore header line
while i<len(pemberton): 
    attributes = pemberton[i].split("\t")
    alternateName = attributes[3]
    repeatStructure = regex.sub('', attributes[9])
    pemberton_dict[alternateName] = repeatStructure
    i+=1

trf_output = open("trf_input.txt.2.30.30.80.10.16.2000.dat")
trf_output = trf_output.readlines()
n = 0 
trf_dict={}
while n< len(trf_output):
    if "Sequence" in trf_output[n]:
        locus = (trf_output[n][10:].rstrip()).split("_")[1]
        n+=7 
        locus_metrics=[]
        while (n<len(trf_output)-1 and trf_output[n]!= '\n'):
            locus_metrics.append(trf_output[n])
            n+=1
        trf_dict[locus]=get_metrics_match(locus_metrics)
    n+=1
with open('filtered_trf_output_debugged.txt', 'w') as convert_file:
     convert_file.write(json.dumps(trf_dict))