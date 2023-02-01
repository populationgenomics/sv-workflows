"""
This script creates a .JSON catalog compatible with ExpansionHunter. The input is a BED file containing loci of pure_repeats, which was outputted from `pure_repeats_catalog.py`

Required packages: json
"""

import json

file= open("intermediate_files/pure_repeats_loci.bed")
file = file.readlines()
eh_catalog =[]

chr9_27573528_27573546_off_target_regions = [
            "chr1:8423620-8424123",
            "chr1:30908133-30908774",
            "chr1:150579390-150580033",
            "chr10:930904-931660",
            "chr10:24952416-24952917",
            "chr10:43304792-43305293",
            "chr10:46336798-46337407",
            "chr10:46816083-46816694",
            "chr10:100999603-101000210",
            "chr11:910565-911217",
            "chr11:2884462-2885264",
            "chr11:44726686-44727485",
            "chr11:67372279-67372872",
            "chr12:459654-460155",
            "chr12:54082197-54082983",
            "chr12:128266836-128267355",
            "chr13:20987998-20988595",
            "chr13:25468389-25469012",
            "chr13:106917207-106917872",
            "chr14:32938886-32939513",
            "chr14:105679947-105680448",
            "chr15:30223043-30223669",
            "chr16:768799-769576",
            "chr16:1614059-1614704",
            "chr16:3134791-3135352",
            "chr16:87602166-87602667",
            "chr17:10198365-10198994",
            "chr17:68205352-68205991",
            "chr17:81891215-81891716",
            "chr18:54223897-54224704",
            "chr19:613299-613854",
            "chr19:878542-879249",
            "chr19:2808202-2808704",
            "chr19:17947549-17948080",
            "chr19:39406514-39407140",
            "chr19:50568087-50568733",
            "chr19:52296834-52297335",
            "chr2:3592013-3592538",
            "chr2:25129074-25129575",
            "chr2:168247028-168247668",
            "chr2:234496868-234497369",
            "chr2:234951698-234952200",
            "chr2:236237107-236237609",
            "chr2:239737824-239738572",
            "chr20:407681-408343",
            "chr20:62064993-62065670",
            "chr21:45286831-45287440",
            "chr22:39994753-39995382",
            "chr22:49117029-49117530",
            "chr22:50627810-50628453",
            "chr3:46845815-46846372",
            "chr3:47578060-47578713",
            "chr3:49554283-49555050",
            "chr3:50124039-50124697",
            "chr3:51706744-51707245",
            "chr3:128497689-128498194",
            "chr3:184379811-184380428",
            "chr4:151408887-151409388",
            "chr5:10761008-10761635",
            "chr5:146458824-146459326",
            "chr5:172454307-172454808",
            "chr5:177370694-177371467",
            "chr6:33633374-33634093",
            "chr6:43021418-43022045",
            "chr6:46490762-46491263",
            "chr6:109483278-109483938",
            "chr8:141308350-141308989",
            "chr9:83817587-83818088",
            "chr9:129568833-129569572",
            "chr9:129665778-129666418",
            "chrX:9785974-9786751",
            "chrX:49879636-49880326",
            "chrX:53625582-53626227",
            "chrX:154341221-154341987",
            "chrX:154390157-154390883"
        ]

for i in file: 
    chr, start, end, motif_length, motif, num_repeats = i.split("\t")
    locus_id = f'{chr}_{start}_{end}'
    locus_structure = f'({motif})*'
    reference_region = f'{chr}:{start}-{end}'
    variant_type ='Repeat'
    entry = {}
    variant_id = locus_id
    entry['LocusId'] = locus_id
    entry['LocusStructure'] = locus_structure
    entry['ReferenceRegion'] = reference_region
    entry['VariantId'] = variant_id
    if reference_region == "chr9:27573528-27573546": #special locus because Illumina catalog specifies Off-Target-Regions, which also changes classification to RareRepeat
        variant_type ='RareRepeat'
    entry['VariantType'] = variant_type 
    if reference_region == "chr9:27573528-27573546": #special locus because Illumina catalog specifies Off-Target-Regions, which also changes classification to RareRepeat
        entry['OfftargetRegions'] = chr9_27573528_27573546_off_target_regions
    eh_catalog.append(entry)

print(len(eh_catalog)) #164847
with open("catalogs/eh_pure_repeats_catalog.json", "w") as eh:
        json.dump(eh_catalog, eh, indent=4)

