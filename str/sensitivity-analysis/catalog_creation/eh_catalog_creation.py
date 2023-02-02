
from str_analysis.convert_gangstr_spec_to_expansion_hunter_variant_catalog import process_variant_catalog

process_variant_catalog("marshfield_regions_gangstr.bed", "marshfield_regions_eh.json")

"""
#Remake gangstr catalog after EH creation to make sure the loci are trimmed (feature of Weisburd's scripts)

from str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec import process_variant_catalog
process_variant_catalog("marshfield_regions_eh.json", "marshfield_regions_gangstr.bed")
"""
