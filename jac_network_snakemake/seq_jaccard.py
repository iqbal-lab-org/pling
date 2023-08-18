import subprocess
from pathlib import Path
from typing import Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from intervaltree import IntervalTree

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

def seq_jaccard(plasmid_1: Path, plasmid_2: Path, prefix: str, plasmid_1_name, plasmid_2_name, identity_threshold=80, length_threshold=200) -> Tuple[str, str]:
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

    for extension in [".1coords", ".1delta", ".delta", ".mcoords", ".mdelta", ".qdiff", ".rdiff", ".report", ".snps"]:
        Path(prefix+extension).unlink()

    coverage_ref = get_coverage(ref_to_block)
    coverage_query = get_coverage(query_to_block)

    if len_ref>len_query:
        jaccard = coverage_query/len_query
    else:
        jaccard = coverage_ref/len_ref

    return jaccard

plasmid_1 = snakemake.params.genome1
plasmid_2 = snakemake.params.genome2
identity_threshold = snakemake.params.identity_threshold
jaccardpath = snakemake.output.jaccard
jaccard = seq_jaccard(snakemake.input.genome_1_fasta, snakemake.input.genome_2_fasta,f"{plasmid_1}~{plasmid_2}", plasmid_1, plasmid_2, identity_threshold)
with open(jaccardpath, 'w+') as f:
    f.write(f"{plasmid_1}\t{plasmid_2}\t{jaccard}\n")
