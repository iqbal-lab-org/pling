import numpy as np
import pandas as pd
import itertools

genomes = snakemake.params.genomes
input = snakemake.input.unimog

sequences = {}
with open(input) as unimog:
    genome = ""
    for line in unimog:
        if line[0]==">":
            genome = line.strip(">\n")
        else:
            sequence=line.strip(")\n").split()
            sequences[genome]=[int(el) for el in sequence]

n = len(genomes)
jaccard_matrix = np.empty((n,n))
jaccard_matrix.fill(np.nan)
with open(gene_jaccard_outfile, "w") as gene_jaccard_fh:
    for genome_1, genome_2 in itertools.combinations(genomes, 2):
        count = 0
        for el in sequences[genome_1]:
            if el in sequences[genome_2]:
                count_1 = sequences[genome_1].count(el)
                count_2 = sequences[genome_2].count(el)
                count = count + min(count_1, count_2)
        jaccard = count/min(len(sequence[genome_1]), len(sequences[genome_2]))
        print(f"{genome_1}\t{genome_2}\t{jaccard}", file=gene_jaccard_fh)
