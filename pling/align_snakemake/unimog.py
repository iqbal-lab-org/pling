from integerise_plasmids import integerise_plasmids
import pandas as pd
import argparse
import os
from pathlib import Path
from pling.utils import read_in_batch_pairs, get_fasta_file_info

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, containment_threshold, identity_threshold, topology_1, topology_2):
    plasmid_1_unimogs, plasmid_2_unimogs, containment_distance, blocks_ref, blocks_query = integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                f"{genome_1}~{genome_2}", genome_1, genome_2, containment_threshold, identity_threshold, topology_1, topology_2)
    unimog = f">{genome_1}~{genome_2}:{genome_1}\n{plasmid_1_unimogs}\n>{genome_1}~{genome_2}:{genome_2}\n{plasmid_2_unimogs}\n"
    return unimog, containment_distance, blocks_ref, blocks_query

def batchwise_unimog(fastafiles, pairs, unimogpath, mappath, containmentpath, identity_threshold, containment_threshold, topologies):
    containments = []
    unimogs = []
    batch_blocks = {}
    for pair in pairs:
        genome_1 = pair[0]
        genome_2 = pair[1]
        genome_1_fasta = fastafiles[genome_1]
        genome_2_fasta = fastafiles[genome_2]
        topology_1 = topologies[genome_1]
        topology_2 = topologies[genome_2]
        unimog, containment, blocks_ref, blocks_query = unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, containment_threshold, identity_threshold, topology_1, topology_2)
        containments.append(f"{genome_1}\t{genome_2}\t{containment}\n")
        if containment<=containment_threshold:
            unimogs.append(unimog)
            blocks = pd.concat([blocks_ref, blocks_query], ignore_index=True)
            batch_blocks[f"{genome_1}~{genome_2}"] = blocks
    with open(unimogpath, 'w') as f:
        for line in unimogs:
            f.write(line)
    with open(mappath, 'w') as f:
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
    with open(containmentpath, 'w') as f:
        for line in containments:
            f.write(line)


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, containment and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--genomes_list", required=True, help="Text file with list of all fasta filepaths")
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--identity_threshold", required=True, type=float, help="Identity threshold for comparison")
    parser.add_argument("--containment_distance", required=True, type=float, help="containment distance threshold value")
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--containment_output", required=True, help="Output path for containment index results")
    parser.add_argument("--unimog_output", required=True, help="Output path for unimog")
    parser.add_argument("--map_output", required=True, help="Output path for map")
    parser.add_argument("--regions", action="store_true", help="Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids.")
    parser.add_argument("--topology", help="File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under \"plasmid\" and one with \"linear\" or \"circular\" as entries under \"topology\". Without this file, pling will asume all plasmids are circular.")

    # Parse the arguments
    args = parser.parse_args()

    fastafiles, fastaext, fastapath = get_fasta_file_info(args.genomes_list)

    pairs=read_in_batch_pairs(f"{args.outputpath}/batches/batch_{args.batch}.txt")
    batch=list({el[0] for el in pairs}.union({el[1] for el in pairs}))

    print(args.regions)
    if args.regions:
        topologies = {plasmid:"region" for plasmid in batch}
    elif args.topology:
        topology_df=pd.read_csv(args.topology, sep="\t")
        topologies = {plasmid:topology_df[topology_df["plasmid"]==plasmid]["topology"].values[0] for plasmid in batch}
    else:
        topologies = {plasmid:"circular" for plasmid in batch}

    batchwise_unimog(fastafiles, pairs, args.unimog_output, args.map_output, args.containment_output, args.identity_threshold, args.containment_distance, topologies)


if __name__ == "__main__":
    main()
