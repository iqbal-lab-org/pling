from integerise_plasmids import integerise_plasmids
from matches import new_integerise_plasmids
import pandas as pd
import argparse
import os
from pathlib import Path

def unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold):
    plasmid_1_unimogs, plasmid_2_unimogs, jaccard_distance, blocks_ref, blocks_query = new_integerise_plasmids(genome_1_fasta, genome_2_fasta,
                                                                f"{genome_1}~{genome_2}", genome_1, genome_2, identity_threshold)
    unimog = f">{genome_1}\n{plasmid_1_unimogs}\n>{genome_2}\n{plasmid_2_unimogs}"
    return unimog, jaccard_distance, blocks_ref, blocks_query

def batchwise_unimog(fastapath, fastaext, pairs, output_dir, jaccardpath, identity_threshold, jaccard_threshold):
    jaccards = []
    for pair in pairs:
        genome_1 = pair[0]
        genome_2 = pair[1]
        genome_1_fasta = f"{fastapath}/{genome_1}{fastaext[genome_1]}"
        genome_2_fasta = f"{fastapath}/{genome_2}{fastaext[genome_2]}"
        unimog, jaccard, blocks_ref, blocks_query = unimogs_to_ilp_core(genome_1_fasta, genome_2_fasta, genome_1, genome_2, identity_threshold)
        jaccards.append(f"{genome_1}\t{genome_2}\t{jaccard}\n")
        if jaccard<=jaccard_threshold:
            filepath = f"{output_dir}/{genome_1}~{genome_2}_align.unimog"
            with open(filepath, 'w+') as f:
                f.write(unimog)
            blocks = pd.concat([blocks_ref, blocks_query])
            blocks.to_csv(f"{output_dir}/{genome_1}~{genome_2}_map.txt", sep="\t", index = False)
    with open(jaccardpath, 'w+') as f:
        for line in jaccards:
            f.write(line)

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, jaccard and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--genomes_list", required=True, help="Text file with list of all fasta filepaths")
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--identity_threshold", required=True, type=float, help="Identity threshold value")
    parser.add_argument("--jaccard_distance", required=True, type=float, help="Jaccard distance threhsold value")
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--jaccard_output", required=True, help="Output path for jaccard")

    # Parse the arguments
    args = parser.parse_args()

    FASTAFILES = [el[0] for el in pd.read_csv(args.genomes_list, header=None).values]
    fastaext = {os.path.splitext(os.path.basename(el))[0]:os.path.splitext(os.path.basename(el))[1] for el in FASTAFILES}
    fastapath = os.path.dirname(FASTAFILES[0])

    pairs=[]
    with open(f"{args.outputpath}/tmp_files/batches/batch_{args.batch}.txt","r") as f:
        for line in f:
            genome1, genome2 = (line.strip("[]\n").split(","))
            pairs.append([genome1[1:-1], genome2[2:-1]])

    output_dir = Path(f"{args.outputpath}/unimogs/batch_{args.batch}")
    output_dir.mkdir(parents=True, exist_ok=True)
    batchwise_unimog(fastapath, fastaext, pairs, output_dir, args.jaccard_output, args.identity_threshold, args.jaccard_distance)


if __name__ == "__main__":
    main()
