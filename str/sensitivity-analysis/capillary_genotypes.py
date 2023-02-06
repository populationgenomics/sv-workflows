import warnings
from itertools import chain, islice
import json


warnings.simplefilter(action='ignore', category=FutureWarning)


def generator_chunks(generator, size):
    """
    Iterates across a generator, returning specifically sized chunks
    Args:
        generator (): any generator or method implementing yield
        size (): size of iterator to return
    Returns:
        a subset of the generator results
    """
    iterator = iter(generator)
    for first in iterator:
        yield list(chain([first], islice(iterator, size - 1)))


file = open("combinedmicrosats_627loci_1048indivs_numRpts.stru.txt")
cpg_hgdp = open("hgdp_cpg_key.json")
sample_dict = json.load(cpg_hgdp)
sample_dict_inverted = {}
for key in sample_dict: 
    sample_dict_inverted[sample_dict[key]] = key
# file = file.readlines()

loci = file.readline().rstrip().split()
# print(len(loci)) #627

motif_length = file.readline().rstrip().split()
# print(len(motif_length)) #627

num_str_regions = file.readline().rstrip().split()
# print(len(num_str_regions)) #627

repeat_regularity = file.readline().rstrip().split()

loci_dict = {
    locus: (motif_len, num_str, regularity)
    for locus, motif_len, num_str, regularity in zip(
        loci, motif_length, num_str_regions, repeat_regularity
    )
}

# print(loci_dict)

with open("capillary_genotypes.csv", 'w', encoding='utf-8') as handle:
    handle.write(
        ','.join(
            [
                'sample_id',
                'cpg_id',
                'locus_id',
                'genotype_1',
                'genotype_2',
                'population_id',
                'population_name',
                'geographic_population',
                'predefined_region',
                'regularity',
                'motif_length',
            ]
        )
        + '\n'
    )

    for line_1, line_2 in generator_chunks(file, 2):

        attributes_1 = line_1.split()
        attributes_2 = line_2.split()
        # checks
        assert attributes_1[:5] == attributes_2[:5]

        if len(str(attributes_1[0])) ==1: 
            sample_id = "HGDP0000"+ attributes_1[0]
        elif len(str(attributes_1[0])) ==2: 
            sample_id = "HGDP000"+attributes_1[0]
        elif len(str(attributes_1[0])) ==3: 
            sample_id = "HGDP00"+attributes_1[0]
        elif len(str(attributes_1[0])) ==4: 
            sample_id = "HGDP0"+attributes_1[0]
        try:
             cpg_id = sample_dict_inverted[sample_id]
        except: 
            print(f'{sample_id} does not map to a CPG id')
            continue
        population_id = attributes_1[1]
        population_name = attributes_1[2]
        geographic_population = attributes_1[3]
        predefined_region = attributes_1[4]

        genotypes_1 = attributes_1[5:]
        genotypes_2 = attributes_2[5:]

        for i, locus in enumerate(loci):

            motif_length, num_str_regions, repeat_regularity = loci_dict[locus]
            handle.write(
                ','.join(
                    [
                        sample_id,
                        cpg_id,
                        locus,
                        genotypes_1[i],
                        genotypes_2[i],
                        population_id,
                        population_name,
                        geographic_population,
                        predefined_region,
                        repeat_regularity,
                        motif_length,
                    ]
                )
                + '\n'
            )
            