from intervaltree import IntervalTree, Interval
from integerise_plasmids import *
from pathlib import Path
import subprocess


class Match:
    def __init__(self, rstart, rend, qstart, qend, strand):
        self.rstart = rstart
        self.rend = rend
        self.qstart = qstart
        self.qend = qend
        self.strand = strand

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __hash__(self):
        return hash((self.rstart, self.rend, self.qstart, self.qend, self.strand))

    def __str__(self):
        return f"({self.rstart}, {self.rend}, {self.qstart}, {self.qend}, {self.strand})"

    def rlen(self):
        return self.rend-self.rstart

    def qlen(self):
        return self.qend-self.qstart

    def indels(self):
        return self.rlen()-self.qlen()

class Matches:
    def __init__(self, list_of_matches): #list_of_matches is a list of Match objects
        self.list = list_of_matches
        self.reference = IntervalTree()
        self.query = IntervalTree()
        for match in self.list:
            self.reference[match.rstart:match.rend] = (match.qstart, match.qend, match.strand)
            self.query[match.qstart:match.qend] = (match.rstart, match.rend, match.strand)

    def __getitem__(self, key):
        return self.list[key]

    def __setitem__(self, key, match):
        old_match = self.list[key]
        self.reference.remove(Interval(old_match.rstart, old_match.rend, (old_match.qstart, old_match.qend, old_match.strand)))
        self.query.remove(Interval(old_match.qstart, old_match.qend, (old_match.rstart, old_match.rend, old_match.strand)))
        if match.rend-match.rstart>0 and match.qend-match.qstart>0:
            self.reference[match.rstart:match.rend] = (match.qstart, match.qend, match.strand)
            self.query[match.qstart:match.qend] = (match.rstart, match.rend, match.strand)
        self.list[key] = match

    def __len__(self):
        return len(self.list)

    def __str__(self):
        matches = [str(match) for match in self.list]
        return '['+', '.join(matches)+']'

    def sort(self, ref_bool): #sort in ascending ref or query start positions, with matches w same start being sorted according to ascending end position
        if ref_bool:
            f_start = lambda match: match.rstart
            f_end = lambda match: match.rend
        else:
            f_start = lambda match: match.qstart
            f_end = lambda match: match.qend
        self.list = sorted(self.list, key=f_start)
        n = len(self.list)
        i=0
        while i < n:
            same_start_matches = [self[i]]
            same_start = True
            j = i + 1
            while same_start and j<n:
                if ref_bool and self[i].rstart==self[j].rstart:
                    same_start_matches.append(self[j])
                elif (not ref_bool) and self[i].qstart==self[j].qstart:
                    same_start_matches.append(self[j])
                else:
                    same_start = False
                j = j+1
            subsort = sorted(same_start_matches, key=f_end)
            for k in range(len(subsort)):
                self.list[i+k] = subsort[k]
            i = i+len(subsort)

    def contain_interval(self, start, end, ref_bool): #find in which matches interval (start, end) is contained in ref/query genome
        matches = []
        if ref_bool:
            interval_matches = list(self.reference[start:end])
            for interval in interval_matches:
                    matches.append(Match(interval.begin, interval.end, interval.data[0], interval.data[1], interval.data[2]))
        else:
            interval_matches = list(self.query[start:end])
            for interval in interval_matches:
                    matches.append(Match(interval.data[0], interval.data[1], interval.begin, interval.end, interval.data[2]))
        return matches

    def split_match(self, match, start, end, ref_bool): #given an interval (start,end) on ref/query genome, split up the match (containing it)
        index = self.list.index(match)
        if ref_bool:
            lhs_len = start - match.rstart
            rhs_len = match.rend - end
            if match.strand == 1:
                lhs_split = Match(match.rstart, start, match.qstart, match.qstart + lhs_len, 1)
                interval = Match(start, end, match.qstart + lhs_len, match.qend - rhs_len, 1)
                rhs_split = Match(end, match.rend, match.qend - rhs_len, match.qend, 1)
            else:
                rhs_split = Match(match.rstart, start, match.qend-lhs_len, match.qend, -1 )
                interval = Match(start, end, match.qstart+rhs_len, match.qend-lhs_len, -1)
                lhs_split = Match(end, match.rend, match.qstart, match.qstart +rhs_len, -1)
        else:
            lhs_len = start - match.qstart
            rhs_len = match.qend - end
            if match.strand == 1:
                lhs_split = Match(match.rstart, match.rstart + lhs_len, match.qstart, start,  1)
                interval = Match(match.rstart + lhs_len, match.rend - rhs_len, start, end, 1)
                rhs_split = Match(match.rend - rhs_len, match.rend, end, match.qend, 1)
            else:
                rhs_split = Match(match.rend-lhs_len, match.rend, match.qstart, start, -1 )
                interval = Match(match.rstart+rhs_len, match.rend-lhs_len, start, end, -1)
                lhs_split = Match(match.rstart, match.rstart +rhs_len, end, match.qend, -1)
        self[index] = lhs_split
        self.list.insert(index+1, interval)
        self.list.insert(index+2, rhs_split)

    def greedy(self, i, roverlap, ref_bool):
        if ref_bool:
            f = lambda x, j: x-self[j].indels() if self[j].indels()>0 else x #if there are indels, overlap might be bigger than the length of the query match
            if self[i].rend-self[i].rstart >= self[i+1].rend-self[i+1].rstart: #lhs match bigger than rhs match, so overlap remains in lhs and is removed from rhs match
                overlap = f(roverlap,i+1)
                if self[i+1].strand == 1:
                    self[i+1]=Match(self[i].rend, self[i+1].rend, self[i+1].qstart+overlap, self[i+1].qend, 1)
                else:
                    self[i+1]=Match(self[i].rend, self[i+1].rend, self[i+1].qstart, self[i+1].qend-overlap, -1)
            else: #rhs match bigger than lhs match, so overlap remains in rhs and is removed from lhs match
                overlap = f(roverlap,i)
                if self[i].qstrand == 1:
                    self[i]=Match(self[i].rstart, self[i+1].rstart, self[i].qstart, self[i].qend-overlap, 1)
                else:
                    self[i]=Match(self[i].rstart, self[i+1].rstart, self[i].qstart+overlap, self[i].qend, -1)
        else:
            f = lambda x, j: x+self[j].indels() if self[j].indels()<0 else x #if there are indels, overlap might be bigger than the length of the query match
            if self[i].qend-self[i].qstart >= self[i+1].qend-self[i+1].qstart: #lhs match bigger than rhs match, so overlap remains in lhs and is removed from rhs match
                overlap = f(roverlap,i+1)
                if self[i+1].strand == 1:
                    self[i+1]=Match(self[i+1].rstart+overlap, self[i+1].rend, self[i].qend, self[i+1].qend, 1)
                else:
                    self[i+1]=Match(self[i+1].rstart, self[i+1].rend-overlap, self[i].qend, self[i+1].qend, -1)
            else: #rhs match bigger than lhs match, so overlap remains in rhs and is removed from lhs match
                overlap = f(roverlap,i)
                if self[i].qstrand == 1:
                    self[i]=Match(self[i].rstart, self[i].rend-overlap, self[i].qstart, self[i+1].qstart, 1)
                else:
                    self[i]=Match(self[i].rstart+overlap, self[i].rend, self[i].qstart, self[i+1].qstart, -1)

    def find_opposite_overlaps(self, i, roverlap, ref_bool):
        overlaps = []
        rstart_1 = self[i].rstart
        rend_1 = self[i].rend
        qstart_1 = self[i].qstart
        qend_1 = self[i].qend
        rstart_2 = self[i+1].rstart
        rend_2 = self[i+1].rend
        qstart_2 = self[i+1].qstart
        qend_2 = self[i+1].qend
        if ref_bool:
            f = lambda x, j: x-self[j].indels() if self[j].indels()>0 else x #if there are indels, overlap might be bigger than the length of the query match
            overlap = f(roverlap,i)
            if self[i].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, qend_1-overlap, qend_1, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, qstart_1, qstart_1+overlap, -1))
            overlap = f(roverlap,i+1)
            if self[i+1].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, qstart_2, qstart_2+overlap, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, qend_2-overlap, qend_2, -1))
        else:
            f = lambda x, j: x+self[j].indels() if self[j].indels()<0 else x #if there are indels, overlap might be bigger than the length of the query match
            overlap = f(roverlap,i)
            if self[i].strand == 1:
                overlaps.append(Match(rend_1-overlap, rend_1, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(rstart_1, rstart_1+overlap, qstart_2, qend_1, -1))
            overlap = f(roverlap,i+1)
            if self[i+1].strand == 1:
                overlaps.append(Match(rstart_2, rstart_2+overlap, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(rend_2-overlap, rend_2, qstart_2, qend_1, -1))
        return overlaps

    def resolve_overlaps(self, overlap_threshold):
        self.sort(True)
        i=0
        finished = False
        while not finished:
            try:
                overlap = self[i].rend-self[i+1].rstart
            except IndexError:
                overlap = 0
                finished = True
            if overlap>overlap_threshold:
                overlap_matches = self.find_opposite_overlaps(i, overlap, True)
                contain_overlap_1 = self.contain_interval(overlap_matches[0].qstart, overlap_matches[0].qend, False)
                contain_overlap_2 = self.contain_interval(overlap_matches[1].qstart, overlap_matches[1].qend, False)
                if len(contain_overlap_1)<3 or len(contain_overlap_2)<3:
                    for match in contain_overlap_1:
                        self.split_match(match, overlap_matches[0].qstart, overlap_matches[0].qend, False)
                    contain_overlap_2 = self.contain_interval(overlap_matches[1].qstart, overlap_matches[1].qend, False)
                    for match in contain_overlap_2:
                        self.split_match(match, overlap_matches[1].qstart, overlap_matches[1].qend, False)
                else:
                    self.greedy(i, overlap, True)
            i = i+1
            #end = self[i].rend

        self.sort(False)
        i=0
        finished = False
        while not finished:
            try:
                overlap = self[i].qend-self[i+1].qstart
            except IndexError:
                overlap = 0
                finished = True
            if overlap>overlap_threshold:
                overlap_matches = self.find_opposite_overlaps(i, overlap, False)
                contain_overlap_1 = self.contain_interval(overlap_matches[0].rstart, overlap_matches[0].rend, True)
                contain_overlap_2 = self.contain_interval(overlap_matches[1].rstart, overlap_matches[1].rend, True)
                if len(contain_overlap_1)<3 or len(contain_overlap_2)<3:
                    for match in contain_overlap_1:
                        self.split_match(match, overlap_matches[0].rstart, overlap_matches[0].rend, True)
                    contain_overlap_2 = self.contain_interval(overlap_matches[1].rstart, overlap_matches[1].rend, True)
                    for match in contain_overlap_2:
                        self.split_match(match, overlap_matches[1].rstart, overlap_matches[1].rend, True)
                else:
                    self.greedy(i, overlap, False)
            i = i+1

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

def new_integerise_plasmids(plasmid_1: Path, plasmid_2: Path, prefix: str, plasmid_1_name, plasmid_2_name, identity_threshold=80, length_threshold=200):
    subprocess.check_call(f"perl -w $(which dnadiff) {plasmid_1} {plasmid_2} -p {prefix} 2>/dev/null", shell=True)
    show_coords_output = subprocess.check_output(f"show-coords -TrcldH -I {identity_threshold} {prefix}.1delta", shell=True).strip().split(b'\n')  # TODO: what about this threshold?

    assert(len(show_coords_output)>0)

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
        og_matches.append(Match(start_ref, end_ref, start_query, end_query, strand_query))

    matches = Matches(og_matches)
    coverage_ref = 0
    coverage_query = 0
    no_matches_available = len_ref==-1
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
        jaccard = coverage_query/len_query
    else:
        jaccard = coverage_ref/len_ref

    return plasmid_1_unimogs, plasmid_2_unimogs, jaccard, blocks_ref, blocks_query

testing = False
if testing == True:
    matches = Matches([Match(1, 1032, 1, 1032, 1), Match(1005, 3385, 550, 2930, 1)])
    overlap_threshold = 200
    length_threshold=0
    matches.resolve_overlaps(overlap_threshold)
    print(matches)
    ref_to_block, query_to_block, max_id = make_interval_tree_w_dups(matches.list, length_threshold)
    populate_interval_tree_with_unmatched_blocks(ref_to_block, 3385, max_id+1, length_threshold)
    populate_interval_tree_with_unmatched_blocks(query_to_block, 2930, len(ref_to_block)+1, length_threshold)
    plasmid_1_unimogs = get_unimog(ref_to_block)
    plasmid_2_unimogs = get_unimog(query_to_block)
    blocks_ref = get_blocks('reference', ref_to_block)
    blocks_query = get_blocks('query', query_to_block)
    print(plasmid_1_unimogs, plasmid_2_unimogs, blocks_ref, blocks_query)
    #plasmid1 = "NZ_LT985234.1"
    #plasmid2 = "NZ_CP062902.1"
    '''
    plasmid2 = "NZ_LR999867.1"
    plasmid1 = "NZ_CP032890.1"
    path = "/home/daria/Documents/projects/INC-plasmids/samples/fastas/incy"
    print(new_integerise_plasmids(f"{path}/{plasmid1}.fna", f"{path}/{plasmid2}.fna", f"{plasmid1}~{plasmid2}", plasmid1, plasmid2, length_threshold=200))'''
