import warnings


warnings.simplefilter(action='ignore', category=FutureWarning)


file = open("combinedmicrosats_627loci_1048indivs_numRpts.stru.txt")
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
                'locus_id',
                'genotype',
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

    for line in file:
        attributes = line.split()
        sample_id = attributes[0]
        population_id = attributes[1]
        population_name = attributes[2]
        geographic_population = attributes[3]
        predefined_region = attributes[4]

        for i, locus in enumerate(loci):

            motif_length, num_str_regions, repeat_regularity = loci_dict[locus]
            handle.write(
                ','.join(
                    [
                        sample_id,
                        locus,
                        attributes[5 + i],
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