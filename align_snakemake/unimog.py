from integerise_plasmids import integerise_plasmids
from matches import new_integerise_plasmids
import pandas as pd

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold):
    plasmid_1_unimogs, plasmid_2_unimogs, jaccard, blocks_ref, blocks_query = new_integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                f"{genome_1}~{genome_2}", genome_1, genome_2, identity_threshold)
    unimog = f">{genome_1}\n{plasmid_1_unimogs}\n>{genome_2}\n{plasmid_2_unimogs}"
    return unimog, jaccard, blocks_ref, blocks_query

identity_threshold = snakemake.params.identity_threshold
unimog, jaccard, blocks_ref, blocks_query = unimogs_to_ilp_core(snakemake.input.genome_1_fasta, snakemake.input.genome_2_fasta, snakemake.params.genome1, snakemake.params.genome2, identity_threshold)
filepath = snakemake.output.unimogs
jaccardpath = snakemake.output.jaccard
with open(filepath, 'w+') as f:
    f.write(unimog)
with open(jaccardpath, 'w+') as f:
    f.write(f"{snakemake.params.genome1}\t{snakemake.params.genome2}\t{jaccard}\n")
blocks = pd.concat([blocks_ref, blocks_query])
blocks.to_csv(snakemake.output.seq_blocks, sep="\t", index = False)
