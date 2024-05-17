"""
call out to the process_variant_catalog function in the str_analysis package
"""

from str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec import (
    process_variant_catalog,
)

process_variant_catalog('catalogs/eh_pure_repeats_catalog.json', 'catalogs/gangstr_pure_repeats_catalog.bed')
