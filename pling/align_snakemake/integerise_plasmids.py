import subprocess
from intervaltree import IntervalTree
from pathlib import Path
from typing import Tuple
import pandas as pd
from matches import *

def make_interval_tree_w_dups(block_coords, length_threshold):
    ref_to_block = IntervalTree()
    query_to_block = IntervalTree()
    max_id = 0
    for i in range(len(block_coords)):
        rstart = block_coords[i].rstart
        rend = block_coords[i].rend
        qstart = block_coords[i].qstart
        qend = block_coords[i].qend
        qstrand = block_coords[i].strand
        if rend-rstart>length_threshold and qend-qstart>length_threshold:
            if len(ref_to_block.envelop(rstart,rend)) != 0:
                query_to_block[qstart:qend] = qstrand * list(ref_to_block[rstart:rend])[0].data
            elif len(query_to_block.envelop(qstart,qend)) != 0:
                ref_to_block[rstart:rend] = qstrand * list(query_to_block[qstart:qend])[0].data
            else:
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

def integerise_plasmids(plasmid_1: Path, plasmid_2: Path, prefix: str, plasmid_1_name, plasmid_2_name, identity_threshold=80, length_threshold=200):
    subprocess.check_call(f"nucmer --diagdiff 20 --breaklen 500  --maxmatch -p {prefix} {plasmid_1} {plasmid_2} && delta-filter -1 {prefix}.delta > {prefix}.1delta", shell=True)
    show_coords_output = subprocess.check_output(f"show-coords -TrcldH -I {identity_threshold} {prefix}.1delta", shell=True).strip().split(b'\n')  # TODO: what about this threshold?
    show_snps_output = subprocess.check_output(f"show-snps -TrH {prefix}.1delta", shell=True).strip().split(b'\n')

    assert(len(show_coords_output)>0)

    indels = []
    if len(show_snps_output)>1:
        for line in show_snps_output:
            split_line = line.split(b'\t')
            rsub = str(split_line[1]).replace("b\'","").replace("\'","")
            qsub = str(split_line[2]).replace("b\'","").replace("\'","")
            rstart = int(split_line[0])
            qstart = int(split_line[3])
            if rsub == ".":
                type = "INS"
                try:
                    extend_indel = indels[-1].rstart==rstart and indels[-1].qstart+indels[-1].len==qstart
                except:
                    extend_indel = False
                if extend_indel:
                    indels[-1].len += 1
                else:
                    indels.append(Indel(rstart,qstart,1,type))
            elif qsub == ".":
                type = "DEL"
                try:
                    extend_indel = indels[-1].qstart==qstart and indels[-1].rstart+indels[-1].len==rstart
                except:
                    extend_indel = False
                if extend_indel:
                    indels[-1].len += 1
                else:
                    indels.append(Indel(rstart,qstart,1,type))

    len_ref = -1
    len_query = -1
    og_matches = []

    for line in show_coords_output:
        split_line = line.split(b'\t')

        try:
            start_ref, end_ref = get_coordinates(split_line, 0)
            start_query, end_query = get_coordinates(split_line, 2)
        except ValueError:
            continue

        len_ref = int(split_line[7])
        len_query = int(split_line[8])
        strand_query = int(split_line[12])
        match_indels = []
        for indel in indels:
            if start_ref<=indel.rstart and end_ref>indel.rend and start_query<=indel.qstart and end_query>indel.qend:
                match_indels.append(indel)
        if end_ref-start_ref>length_threshold and end_query-start_query>length_threshold:
            og_matches.append(Match(start_ref, end_ref, start_query, end_query, strand_query, match_indels))

    matches = Matches(og_matches)
    coverage_ref = 0
    coverage_query = 0
    no_matches_available = len(og_matches)==0
    if no_matches_available:
        plasmid_1_unimogs = "1 )"
        plasmid_2_unimogs = "2 )"
        blocks_ref = pd.DataFrame({"Plasmid":[], "Block_ID":[], "Start":[], "End":[]})
        blocks_query = pd.DataFrame({"Plasmid":[], "Block_ID":[], "Start":[], "End":[]})
    else:
        overlap_threshold = 0
        coverage_ref = get_coverage(matches.reference)
        coverage_query = get_coverage(matches.query)
        matches.resolve_overlaps(overlap_threshold)
        ref_to_block, query_to_block, max_id = make_interval_tree_w_dups(matches.list, length_threshold)
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
        jaccard_similarity = coverage_query/len_query
    else:
        jaccard_similarity = coverage_ref/len_ref
    jaccard_distance = 1-jaccard_similarity

    return plasmid_1_unimogs, plasmid_2_unimogs, jaccard_distance, blocks_ref, blocks_query

testing = False

if testing == True:
    #plasmid1 = "NZ_LT985234.1"
    #plasmid2 = "NZ_CP062902.1"
    #plasmid1 = "NZ_LR999867.1"
    #plasmid2 = "NZ_CP032890.1"
    plasmid1 = "cpe041_12"
    plasmid2 = "cpe024_2"
    #path = "/home/daria/Documents/projects/INC-plasmids/samples/fastas/incy"
    path = "/home/daria/Documents/projects/addenbrookes/same_assembler"
    print(integerise_plasmids(f"{path}/{plasmid1}.fna", f"{path}/{plasmid2}.fna", f"{plasmid1}~{plasmid2}", plasmid1, plasmid2))
