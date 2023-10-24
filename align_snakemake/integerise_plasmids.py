import subprocess
from intervaltree import IntervalTree
from pathlib import Path
from typing import Tuple
import pandas as pd

def split_overlaps_naive(matches, r, q, length_threshold):
    data = []
    rstart = f"{r}start"
    rend = f"{r}end"
    qstart = f"{q}start"
    qend = f"{q}end"
    overlap_bool = False
    for i in range(len(matches)-1):
        overlap = matches[i][rend]-matches[i+1][rstart]
        if overlap>length_threshold:
            if matches[i]["qstrand"] == 1:
                data.append({rstart:matches[i][rstart], rend:matches[i+1][rstart], qstart:matches[i][qstart], qend:matches[i][qend]-overlap, "qstrand":matches[i]["qstrand"]})
                data.append({rstart:matches[i+1][rstart], rend:matches[i][rend], qstart:matches[i][qend]-overlap, qend:matches[i][qend], "qstrand":matches[i]["qstrand"]})
            else:
                data.append({rstart:matches[i][rstart], rend:matches[i+1][rstart], qstart:matches[i][qstart]+overlap, qend:matches[i][qend], "qstrand":matches[i]["qstrand"]})
                data.append({rstart:matches[i+1][rstart], rend:matches[i][rend], qstart:matches[i][qstart], qend:matches[i][qstart]+overlap, "qstrand":matches[i]["qstrand"]})
            if matches[i+1]["qstrand"] == 1:
                data.append({rstart:matches[i+1][rstart], rend:matches[i][rend], qstart:matches[i+1][qstart], qend:matches[i+1][qstart]+overlap, "qstrand":matches[i+1]["qstrand"]})
                data.append({rstart:matches[i][rend], rend:matches[i+1][rend], qstart:matches[i+1][qstart]+overlap, qend:matches[i+1][qend], "qstrand":matches[i+1]["qstrand"]})
            else:
                data.append({rstart:matches[i+1][rstart], rend:matches[i][rend], qstart:matches[i+1][qend]-overlap, qend:matches[i+1][qend], "qstrand":matches[i+1]["qstrand"]})
                data.append({rstart:matches[i][rend], rend:matches[i+1][rend], qstart:matches[i+1][qstart], qend:matches[i+1][qend]-overlap, "qstrand":matches[i+1]["qstrand"]})
        elif i == len(matches)-2:
            data.append({rstart:matches[i][rstart], rend:matches[i][rend], qstart:matches[i][qstart], qend:matches[i][qend], "qstrand":matches[i]["qstrand"]})
            data.append({rstart:matches[i+1][rstart], rend:matches[i+1][rend], qstart:matches[i+1][qstart], qend:matches[i+1][qend], "qstrand":matches[i+1]["qstrand"]})
        else:
            data.append({rstart:matches[i][rstart], rend:matches[i][rend], qstart:matches[i][qstart], qend:matches[i][qend], "qstrand":matches[i]["qstrand"]})
    return data

def split_overlaps_greedy(og_matches, r, q, length_threshold):
    matches = [match for match in og_matches]
    rstart = f"{r}start"
    rend = f"{r}end"
    qstart = f"{q}start"
    qend = f"{q}end"
    overlap_bool = False
    for i in range(len(matches)-1):
        overlap = matches[i][rend]-matches[i+1][rstart]
        if overlap>length_threshold:
            if matches[i][rend]-matches[i][rstart] >= matches[i+1][rend]-matches[i+1][rstart]: #lhs match bigger than rhs match, so overlap remains in lhs and is removed from rhs match
                if matches[i+1]["qstrand"] == 1:
                    matches[i+1]={rstart:matches[i][rend], rend:matches[i+1][rend], qstart:matches[i+1][qstart]+overlap, qend:matches[i+1][qend], "qstrand":matches[i+1]["qstrand"]}
                else:
                    matches[i+1]={rstart:matches[i][rend], rend:matches[i+1][rend], qstart:matches[i+1][qstart], qend:matches[i+1][qend]-overlap, "qstrand":matches[i+1]["qstrand"]}
            else: #rhs match bigger than lhs match, so overlap remains in rhs and is removed from lhs match
                if matches[i]["qstrand"] == 1:
                    matches[i]={rstart:matches[i][rstart], rend:matches[i+1][rstart], qstart:matches[i][qstart], qend:matches[i][qend]-overlap, "qstrand":matches[i]["qstrand"]}
                else:
                    matches[i]={rstart:matches[i][rstart], rend:matches[i+1][rstart], qstart:matches[i][qstart]+overlap, qend:matches[i][qend], "qstrand":matches[i]["qstrand"]}
    return matches

def make_interval_tree_w_dups(block_coords):
    ref_to_block = IntervalTree()
    query_to_block = IntervalTree()
    max_id = 0
    for i in range(len(block_coords)):
        rstart = block_coords[i]["rstart"]
        rend = block_coords[i]["rend"]
        qstart = block_coords[i]["qstart"]
        qend = block_coords[i]["qend"]
        qstrand = block_coords[i]["qstrand"]
        if len(ref_to_block.envelop(rstart,rend)) != 0:
            query_to_block[qstart:qend] = qstrand * list(ref_to_block[rstart:rend])[0].data
        elif len(query_to_block.envelop(qstart,qend)) != 0:
            ref_to_block[rstart:rend] = qstrand * list(query_to_block[qstart:qend])[0].data
        else:
            max_id=max_id+1
            ref_to_block[rstart:rend] = max_id
            query_to_block[qstart:qend] = qstrand * max_id
    return ref_to_block, query_to_block, max_id

def make_interval_tree(block_coords):
    ref_to_block = IntervalTree()
    query_to_block = IntervalTree()
    max_id = 0
    for i in range(len(block_coords)):
        rstart = block_coords[i]["rstart"]
        rend = block_coords[i]["rend"]
        qstart = block_coords[i]["qstart"]
        qend = block_coords[i]["qend"]
        qstrand = block_coords[i]["qstrand"]
        max_id=max_id+1
        ref_to_block[rstart:rend] = max_id
        query_to_block[qstart:qend] = qstrand * max_id
    return ref_to_block, query_to_block, max_id

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

def get_coverage(tree):
    merge_tree = IntervalTree(list(tree))
    merge_tree.merge_overlaps()
    coverage = 0
    for interval in merge_tree:
        coverage = coverage + (interval.end - interval.begin)
    return coverage

def integerise_plasmids(plasmid_1: Path, plasmid_2: Path, prefix: str, plasmid_1_name, plasmid_2_name, identity_threshold=80, length_threshold=200) -> Tuple[str, str]:
    subprocess.check_call(f"perl -w $(which dnadiff) {plasmid_1} {plasmid_2} -p {prefix} 2>/dev/null", shell=True)
    show_coords_output = subprocess.check_output(f"show-coords -TrcldH -I {identity_threshold} {prefix}.1delta", shell=True).strip().split(b'\n')  # TODO: what about this threshold?

    assert(len(show_coords_output)>0)

    len_ref = -1
    len_query = -1
    og_matches = []
    ref_coverage = IntervalTree()
    query_coverage = IntervalTree()
    for line in show_coords_output:
        split_line = line.split(b'\t')

        try:
            start_ref, end_ref = get_coordinates(split_line, 0)
            start_query, end_query = get_coordinates(split_line, 2)
        except ValueError:
            continue
        ref_coverage[start_ref:end_ref] = 0
        query_coverage[start_query:end_query] = 0
        len_ref = int(split_line[7])
        len_query = int(split_line[8])
        strand_query = int(split_line[12])
        og_matches.append({"rstart":start_ref, "rend":end_ref, "qstart":start_query, "qend":end_query, "qstrand":strand_query})


    coverage_ref = 0
    coverage_query = 0
    no_matches_available = len_ref==-1
    if no_matches_available:
        plasmid_1_unimogs = "1 )"
        plasmid_2_unimogs = "2 )"
        blocks_ref = pd.DataFrame({"Plasmid":[], "Block_ID":[], "Start":[], "End":[]})
        blocks_query = pd.DataFrame({"Plasmid":[], "Block_ID":[], "Start":[], "End":[]})
    else:
        coverage_ref = get_coverage(ref_coverage)
        coverage_query = get_coverage(query_coverage)
        ref_split = split_overlaps_greedy(og_matches, "r", "q", length_threshold)
        blub = sorted(ref_split, key=lambda match: match["qstart"])
        block_coords = split_overlaps_greedy(blub, "q", "r", length_threshold)
        ref_to_block, query_to_block, max_id = make_interval_tree(block_coords)
        populate_interval_tree_with_unmatched_blocks(ref_to_block, len_ref, max_id+1, length_threshold)
        populate_interval_tree_with_unmatched_blocks(query_to_block, len_query, len(ref_to_block)+1, length_threshold)
        plasmid_1_unimogs = get_unimog(ref_to_block)
        plasmid_2_unimogs = get_unimog(query_to_block)
        blocks_ref = get_blocks(plasmid_1_name, ref_to_block)
        blocks_query = get_blocks(plasmid_2_name, query_to_block)

    for extension in [".1coords", ".1delta", ".delta", ".mcoords", ".mdelta", ".qdiff", ".rdiff", ".report", ".snps", ".unqry", ".unref"]:
        try:
            Path(prefix+extension).unlink()
        except:
            pass

    if len_ref>len_query:
        jaccard = coverage_query/len_query
    else:
        jaccard = coverage_ref/len_ref

    return plasmid_1_unimogs, plasmid_2_unimogs, jaccard, blocks_ref, blocks_query

testing = False

if testing == True:
    matches = [{"rstart":1, "rend":1032, "qstart":1, "qend":1032, "qstrand":1}, {"rstart":1005, "rend":3385, "qstart":550, "qend":2930, "qstrand":1}]

    ref_split = split_overlaps_greedy(matches, "r", "q", 200)
    sorting = lambda match: match["qstart"]
    blub = sorted(ref_split, key=sorting)
    print(blub)
    block_coords = split_overlaps_greedy(blub, "q", "r", 200)
    print(block_coords)
    ref_to_block, query_to_block, max_id = make_interval_tree(block_coords)
    populate_interval_tree_with_unmatched_blocks(ref_to_block, 3385, max_id+1, 200)
    populate_interval_tree_with_unmatched_blocks(query_to_block, 2930, len(ref_to_block)+1, 200)
    #print(sorted(ref_to_block))
    #print(sorted(query_to_block))
    print(get_unimog(ref_to_block))
    print(get_unimog(query_to_block))
