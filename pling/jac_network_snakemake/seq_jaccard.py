import subprocess
from pathlib import Path
from typing import Tuple
from intervaltree import IntervalTree
import argparse
from pling.utils import read_in_batch_pairs, get_fasta_file_info

def get_coverage(tree):
    tree.merge_overlaps()
    coverage = 0
    for interval in tree:
        coverage = coverage + (interval.end - interval.begin)
    return coverage

def get_coordinates(split_line, index):
    start = int(split_line[index])
    end = int(split_line[index+1])
    if start > end:
        start, end = end, start
    end+=1
    return start, end

def get_sequence_containment_distance(plasmid_1: Path, plasmid_2: Path, prefix: str, identity_threshold=80) -> Tuple[str, str]:
    subprocess.check_call(f"perl -w $(which dnadiff) {plasmid_1} {plasmid_2} -p {prefix} 2>/dev/null", shell=True)
    show_coords_output = subprocess.check_output(f"show-coords -TrcldH -I {identity_threshold} {prefix}.1delta", shell=True).strip().split(b'\n')  # TODO: what about this threshold?

    assert(len(show_coords_output)>0)

    len_ref = -1
    len_query = -1
    ref_to_block = IntervalTree()
    query_to_block = IntervalTree()
    for block, line in enumerate(show_coords_output):
        split_line = line.split(b'\t')

        try:
            start_ref, end_ref = get_coordinates(split_line, 0)
            start_query, end_query = get_coordinates(split_line, 2)
        except ValueError:
            continue

        len_ref = int(split_line[7])
        len_query = int(split_line[8])
        ref_to_block[start_ref:end_ref] = block
        query_to_block[start_query:end_query] = block

    for extension in [".1coords", ".1delta", ".delta", ".mcoords", ".mdelta", ".qdiff", ".rdiff", ".report", ".snps", ".unqry", ".unref"]:
        try:
            Path(prefix+extension).unlink()
        except:
            pass

    coverage_ref = get_coverage(ref_to_block)
    coverage_query = get_coverage(query_to_block)

    if len_ref>len_query:
        containment_similarity = coverage_query/len_query
    else:
        containment_similarity = coverage_ref/len_ref

    containment_distance = 1-containment_similarity
    return containment_distance

def batchwise_containment(fastafiles, pairs, containmentpath, identity_threshold):
    containments = []
    for pair in pairs:
        genome_1 = pair[0]
        genome_2 = pair[1]
        genome_1_fasta = fastafiles[genome_1]
        genome_2_fasta = fastafiles[genome_2]
        containment_distance = get_sequence_containment_distance(genome_1_fasta, genome_2_fasta, f"{genome_1}~{genome_2}", identity_threshold)
        containments.append(f"{genome_1}\t{genome_2}\t{containment_distance}\n")
    with open(containmentpath, 'w') as f:
        for line in containments:
            f.write(line)


def main(args):
    fastafiles, fastaext, fastapath = get_fasta_file_info(args.genomes_list)

    pairs=read_in_batch_pairs(f"{args.outputpath}/batches/batch_{args.batch}.txt")

    batchwise_containment(fastafiles, pairs, args.containment_output, args.identity_threshold)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate containment index for genome sequences.')

    parser.add_argument("--genomes_list", required=True, help="Text file with list of all fasta filepaths")
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--identity_threshold", required=True, type=float, help="Identity threshold for comparison")
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--containment_output", required=True, help="Output path for containment index results")

    args = parser.parse_args()
    main(args)
