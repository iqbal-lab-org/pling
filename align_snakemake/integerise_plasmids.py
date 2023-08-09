import subprocess
from intervaltree import IntervalTree
from pathlib import Path
from typing import Tuple
import pandas as pd

def populate_interval_tree_with_unmatched_blocks(interval_tree, total_length, block_index, length_threshold):
    pos = 1
    while pos <= total_length:
        no_matches = len(interval_tree[pos])==0
        if no_matches:
            start_pos = pos
            while len(interval_tree[pos])==0 and pos<=total_length:
                pos+=1
            end_pos = pos

            too_short_unmatched_block = end_pos-start_pos <= length_threshold
            if not too_short_unmatched_block:
                interval_tree[start_pos:end_pos]=block_index
                block_index+=1
        else:
            pos+=1


def get_unimog(interval_tree):
    intervals=[]
    for interval in sorted(interval_tree):
        intervals.append(str(interval.data))
    intervals.append(")")
    return ' '.join(intervals)

def get_blocks(plasmid, interval_tree):
    data = {"Plasmid":[], "Block_ID":[], "Start":[], "End":[]}
    for interval in sorted(interval_tree):
        data["Plasmid"].append(plasmid)
        data["Block_ID"].append(interval.data)
        data["Start"].append(interval.begin)
        data["End"].append(interval.end)
    blocks = pd.DataFrame(data)
    return blocks

def get_coordinates(split_line, index):
    start = int(split_line[index])
    end = int(split_line[index+1])
    if start > end:
        start, end = end, start
    end+=1
    return start, end


def integerise_plasmids(plasmid_1: Path, plasmid_2: Path, prefix: str, plasmid_1_name, plasmid_2_name, identity_threshold=80, length_threshold=200) -> Tuple[str, str]:
    subprocess.check_call(f"perl -w $(which dnadiff) {plasmid_1} {plasmid_2} -p {prefix} 2>/dev/null", shell=True)
    show_coords_output = subprocess.check_output(f"show-coords -TrcldH -I {identity_threshold} {prefix}.1delta", shell=True).strip().split(b'\n')  # TODO: what about this threshold?

    assert(len(show_coords_output)>0)

    len_ref = -1
    len_query = -1
    ref_to_block = IntervalTree()
    query_to_block = IntervalTree()
    data = {"Plasmids":[], "Block_ID":[], "Start":[], "Stop":[]}
    coverage_ref = 0
    coverage_query = 0
    for block, line in enumerate(show_coords_output):
        split_line = line.split(b'\t')
        block=block+1  # avoids using block 0, which can't be represented as -0 and 0

        try:
            start_ref, end_ref = get_coordinates(split_line, 0)
            start_query, end_query = get_coordinates(split_line, 2)
        except ValueError:
            continue

        len_ref = int(split_line[7])
        len_query = int(split_line[8])
        strand_ref = int(split_line[11])
        strand_query = int(split_line[12])
        coverage_ref = coverage_ref + (end_ref-start_ref)
        coverage_query = coverage_query + (end_query-start_query)
        ref_to_block[start_ref:end_ref] = block*strand_ref
        query_to_block[start_query:end_query] = block*strand_query

    no_matches_available = len_ref==-1
    if no_matches_available:
        plasmid_1_unimogs = "1 )"
        plasmid_2_unimogs = "2 )"
    else:
        populate_interval_tree_with_unmatched_blocks(ref_to_block, len_ref, len(ref_to_block)+len(query_to_block)+1, length_threshold)
        populate_interval_tree_with_unmatched_blocks(query_to_block, len_query, len(ref_to_block)+len(query_to_block)+1, length_threshold)
        plasmid_1_unimogs = get_unimog(ref_to_block)
        plasmid_2_unimogs = get_unimog(query_to_block)
        blocks_ref = get_blocks(plasmid_1_name, ref_to_block)
        blocks_query = get_blocks(plasmid_2_name, query_to_block)

    for extension in [".1coords", ".1delta", ".delta", ".mcoords", ".mdelta", ".qdiff", ".rdiff", ".report", ".snps"]:
        Path(prefix+extension).unlink()

    if len_ref>len_query:
        jaccard = coverage_query/len_query
    else:
        jaccard = coverage_ref/len_ref

    return plasmid_1_unimogs, plasmid_2_unimogs, jaccard, blocks_ref, blocks_query
