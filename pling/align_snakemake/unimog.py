from integerise_plasmids import integerise_plasmids
from matches import new_integerise_plasmids
import pandas as pd
import argparse

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold):
    plasmid_1_unimogs, plasmid_2_unimogs, jaccard_distance, blocks_ref, blocks_query = new_integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                f"{genome_1}~{genome_2}", genome_1, genome_2, identity_threshold)
    unimog = f">{genome_1}\n{plasmid_1_unimogs}\n>{genome_2}\n{plasmid_2_unimogs}"
    return unimog, jaccard_distance, blocks_ref, blocks_query


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, jaccard and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--genome_1_fasta", required=True, help="Path to genome 1 fasta file")
    parser.add_argument("--genome_2_fasta", required=True, help="Path to genome 2 fasta file")
    parser.add_argument("--genome1", required=True, help="Parameter for genome 1")
    parser.add_argument("--genome2", required=True, help="Parameter for genome 2")
    parser.add_argument("--identity_threshold", required=True, type=float, help="Identity threshold value")
    parser.add_argument("--unimogs_output", required=True, help="Output path for unimogs")
    parser.add_argument("--jaccard_output", required=True, help="Output path for jaccard")
    parser.add_argument("--seq_blocks_output", required=True, help="Output path for sequence blocks")

    # Parse the arguments
    args = parser.parse_args()

    # Use the arguments
    unimog, jaccard_distance, blocks_ref, blocks_query = unimogs_to_ilp_core(args.genome_1_fasta, args.genome_2_fasta,
                                                                    args.genome1, args.genome2, args.identity_threshold)
    with open(args.unimogs_output, 'w') as f:
        f.write(unimog)

    with open(args.jaccard_output, 'w') as f:
        f.write(f"{args.genome1}\t{args.genome2}\t{jaccard_distance}\n")

    blocks = pd.concat([blocks_ref, blocks_query])
    blocks.to_csv(args.seq_blocks_output, sep="\t", index=False)


if __name__ == "__main__":
    main()

