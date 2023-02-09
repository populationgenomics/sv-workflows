from convert_gangstr_spec_to_expansion_hunter_variant_catalog_modified import process_variant_catalog

process_variant_catalog("catalogs/marshfield_regions_gangstr_untrimmed_coordinates.bed", "catalogs/marshfield_regions_eh_trimmed_coordinates.json")



#Remake gangstr catalog after EH creation to make sure the loci are trimmed (feature of Weisburd's scripts)

from convert_expansion_hunter_variant_catalog_to_gangstr_spec_modified import process_variant_catalog
process_variant_catalog("catalogs/marshfield_regions_eh_trimmed_coordinates.json", "catalogs/marshfield_regions_gangstr_trimmed_coordinates.bed")


