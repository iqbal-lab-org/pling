from integerise_plasmids import integerise_plasmids
import pandas as pd
import argparse
import os
from pathlib import Path
from pling.utils import read_in_batch_pairs, get_fasta_file_info

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold):
    plasmid_1_unimogs, plasmid_2_unimogs, jaccard_distance, blocks_ref, blocks_query = integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                f"{genome_1}~{genome_2}", genome_1, genome_2, identity_threshold)
    unimog = f">{genome_1}~{genome_2}:{genome_1}\n{plasmid_1_unimogs}\n>{genome_1}~{genome_2}:{genome_2}\n{plasmid_2_unimogs}\n"
    return unimog, jaccard_distance, blocks_ref, blocks_query

def batchwise_unimog(fastapath, fastaext, pairs, unimogpath, mappath, jaccardpath, identity_threshold, jaccard_threshold):
    jaccards = []
    unimogs = []
    batch_blocks = {}
    for pair in pairs:
        genome_1 = pair[0]
        genome_2 = pair[1]
        genome_1_fasta = f"{fastapath}/{genome_1}{fastaext[genome_1]}"
        genome_2_fasta = f"{fastapath}/{genome_2}{fastaext[genome_2]}"
        unimog, jaccard, blocks_ref, blocks_query = unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold)
        jaccards.append(f"{genome_1}\t{genome_2}\t{jaccard}\n")
        if jaccard<=jaccard_threshold:
            unimogs.append(unimog)
            blocks = pd.concat([blocks_ref, blocks_query], ignore_index=True)
            batch_blocks[f"{genome_1}~{genome_2}"] = blocks
    with open(unimogpath, 'w+') as f:
        for line in unimogs:
            f.write(line)
    with open(mappath, 'w+') as f:
        f.write("\tPlasmid\tBlock_ID\tStart\tEnd\n")
        for key in batch_blocks.keys():
            f.write(f"{key}")
            blocks = batch_blocks[key]
            for index in blocks.index:
                plasmid = blocks.loc[index, "Plasmid"]
                block_id = blocks.loc[index, "Block_ID"]
                start = blocks.loc[index, "Start"]
                end = blocks.loc[index, "End"]
                f.write(f"\t{plasmid}\t{block_id}\t{start}\t{end}\n")
    with open(jaccardpath, 'w+') as f:
        for line in jaccards:
            f.write(line)


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, jaccard and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--genomes_list", required=True, help="Text file with list of all fasta filepaths")
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--identity_threshold", required=True, type=float, help="Identity threshold for comparison")
    parser.add_argument("--jaccard_distance", required=True, type=float, help="Jaccard distance threshold value")
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--jaccard_output", required=True, help="Output path for Jaccard index results")
    parser.add_argument("--unimog_output", required=True, help="Output path for unimog")
    parser.add_argument("--map_output", required=True, help="Output path for map")

    # Parse the arguments
    args = parser.parse_args()

    fastafiles, fastaext, fastapath = get_fasta_file_info(args.genomes_list)

    pairs=read_in_batch_pairs(f"{args.outputpath}/batches/batch_{args.batch}.txt")

    batchwise_unimog(fastapath, fastaext, pairs, args.unimog_output, args.map_output, args.jaccard_output, args.identity_threshold, args.jaccard_distance)


if __name__ == "__main__":
    main()
