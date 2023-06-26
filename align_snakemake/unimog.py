from integerise_plasmids import integerise_plasmids

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2):
    plasmid_1_unimogs, plasmid_2_unimogs = integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                prefix=f"{genome_1}~{genome_2}")
    genome_1_pair_identifier = f"{genome_1}_from_{genome_1}_vs_{genome_2}"
    genome_2_pair_identifier = f"{genome_2}_from_{genome_1}_vs_{genome_2}"
    unimog = f">{genome_1_pair_identifier}\n{plasmid_1_unimogs}\n>{genome_2_pair_identifier}\n{plasmid_2_unimogs}"
    return unimog

unimog = unimogs_to_ilp_core(snakemake.input.genome_1_fasta, snakemake.input.genome_2_fasta, snakemake.params.genome1, snakemake.params.genome2)
filepath = snakemake.output.unimogs
with open(filepath, 'w+') as f:
    f.write(unimog)
