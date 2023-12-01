from intervaltree import IntervalTree, Interval
from integerise_plasmids import *
from pathlib import Path
import subprocess


class Indel:
    def __init__(self, rstart, qstart, len, type):
        self.type = type
        self.len = len
        self.rstart = rstart
        self.qstart = qstart
        self.rend = rstart+len
        self.qend = qstart+len

    def __str__(self):
        return f"({self.rstart}, {self.rend}, {self.qstart}, {self.qend}, {self.type})"

class Match:
    def __init__(self, rstart, rend, qstart, qend, strand, indels=[]):
        self.rstart = rstart
        self.rend = rend
        self.qstart = qstart
        self.qend = qend
        self.strand = strand
        self.indels = indels

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

    def projection(self, coord, ref_bool): #if ref_bool=True, project from reference to query, else query to reference
        if ref_bool:
            dist = coord - self.rstart
            for indel in self.indels:
                if indel.rstart<=coord<indel.rend:
                    return indel.qstart
                elif indel.rstart<coord:
                    if indel.type == "INS":
                        dist += indel.len
                    elif indel.type == "DEL":
                        dist -= indel.len
            if self.strand == 1:
                projected_coord = self.qstart + dist
            else:
                projected_coord = self.qend - dist
        else:
            dist = coord - self.qstart
            for indel in self.indels:
                if indel.qstart<=coord<indel.qend:
                    return indel.rstart
                elif indel.qstart<coord:
                    if indel.type == "INS":
                        dist += indel.len
                    elif indel.type == "DEL":
                        dist -= indel.len
            if self.strand == 1:
                projected_coord = self.rstart + dist
            else:
                projected_coord = self.rend - dist
        return projected_coord

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

    def insert(self, key, match):
        self.list.insert(key, match)
        if match.rend-match.rstart>0 and match.qend-match.qstart>0:
            self.reference[match.rstart:match.rend]= (match.qstart, match.qend, match.strand)
            self.query[match.qstart:match.qend] = (match.rstart, match.rend, match.strand)

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

    def purge_null_intervals(self):
        purged = [el for el in self.list if el.rstart!=el.rend and el.qstart!=el.qend]
        self.list = purged

    def contain_interval(self, start, end, ref_bool): #find in which matches interval (start, end) is contained in ref/query genome
        matches = []
        if ref_bool:
            interval_matches = list(self.reference[start:end])
            for interval in interval_matches:
                if interval.begin <=start and interval.end >=end:
                    matches.append(Match(interval.begin, interval.end, interval.data[0], interval.data[1], interval.data[2]))
        else:
            interval_matches = list(self.query[start:end])
            for interval in interval_matches:
                if interval.begin <=start and interval.end >=end:
                    matches.append(Match(interval.data[0], interval.data[1], interval.begin, interval.end, interval.data[2]))
        return matches

    def split_match(self, match, start, end, ref_bool): #given an interval (start,end) on ref/query genome, split up the match (containing it)
        index = self.list.index(match)
        projected_start = match.projection(start,ref_bool)
        projected_end = match.projection(end,ref_bool)
        if ref_bool:
            if match.strand == 1:
                lhs_split = Match(match.rstart, start, match.qstart, projected_start, 1)
                interval = Match(start, end, projected_start, projected_end, 1)
                rhs_split = Match(end, match.rend, projected_end, match.qend, 1)
            else:
                rhs_split = Match(match.rstart, start, projected_start, match.qend, -1 )
                interval = Match(start, end, projected_end, projected_start, -1)
                lhs_split = Match(end, match.rend, match.qstart, projected_end, -1)
        else:
            if match.strand == 1:
                lhs_split = Match(match.rstart, projected_start, match.qstart, start,  1)
                interval = Match(projected_start, projected_end, start, end, 1)
                rhs_split = Match(projected_end, match.rend, end, match.qend, 1)
            else:
                rhs_split = Match(projected_start, match.rend, match.qstart, start, -1 )
                interval = Match(projected_end, projected_start, start, end, -1)
                lhs_split = Match(match.rstart, projected_end, end, match.qend, -1)
        if lhs_split.rend != lhs_split.rstart or lhs_split.qend!=lhs_split.qstart:
            self[index] = lhs_split #add even if null interval bc otherwise will break the walk from left to right -- null interval from here will be removed later
            self.insert(index+1, interval)
            if rhs_split.rend != rhs_split.rstart or rhs_split.qend!=rhs_split.qstart: #don't add null interval
                self.insert(index+2, rhs_split)
        else:
            self[index] = interval #add even if null interval bc otherwise will break the walk from left to right -- null interval from here will be removed later
            if rhs_split.rend != rhs_split.rstart or rhs_split.qend!=rhs_split.qstart: #don't add null interval
                self.insert(index+1, rhs_split)

    def find_opposite_overlaps(self, i, ref_bool):
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
            if rend_1>rend_2: #check for containment
                rend_1 = self[i+1].rend
            projected_rstart_2 = self[i].projection(rstart_2, ref_bool)
            projected_rend_1 = self[i].projection(rend_1, ref_bool)
            if self[i].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, projected_rstart_2, projected_rend_1, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, projected_rend_1, projected_rstart_2, -1))
            projected_rend_1 = self[i+1].projection(rend_1,ref_bool)
            projected_rstart_2 = self[i+1].projection(rstart_2, ref_bool)
            if self[i+1].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, projected_rstart_2, projected_rend_1, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, projected_rend_1, projected_rstart_2, -1))
        else:
            if qend_1>qend_2:
                qend_1 = self[i+1].qend
            projected_qstart_2 = self[i].projection(qstart_2, ref_bool)
            projected_qend_1 = self[i].projection(qend_1, ref_bool)
            if self[i].strand == 1:
                overlaps.append(Match(projected_qstart_2, projected_qend_1, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(projected_qend_1, projected_qstart_2, qstart_2, qend_1, -1))
            projected_qend_1 = self[i+1].projection(qend_1,ref_bool)
            projected_qstart_2 = self[i+1].projection(qstart_2, ref_bool)
            if self[i+1].strand == 1:
                overlaps.append(Match(projected_qstart_2, projected_qend_1, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(projected_qend_1, projected_qstart_2, qstart_2, qend_1, -1))
        return overlaps

    def resolve_overlaps(self, overlap_threshold):
        self.sort(False)
        i=0
        finished = False
        while not finished:
            overlap = 0
            try:
                overlap = self[i].qend-self[i+1].qstart
                not_null = self[i].qstart!=self[i].qend and self[i+1].qstart!=self[i+1].qend #ignore null intervals
                containment = self[i+1].qend<self[i].qend
            except IndexError:
                finished = True
                not_null = False
                containment = False
            if not_null and overlap>overlap_threshold:
                overlap_matches = self.find_opposite_overlaps(i, False)
                contain_overlap_1 = self.contain_interval(overlap_matches[0].rstart, overlap_matches[0].rend, True)
                for match in contain_overlap_1:
                    self.split_match(match, overlap_matches[0].rstart, overlap_matches[0].rend, True)
                contain_overlap_2 = self.contain_interval(overlap_matches[1].rstart, overlap_matches[1].rend, True)
                for match in contain_overlap_2:
                    self.split_match(match, overlap_matches[1].rstart, overlap_matches[1].rend, True)
                if containment:
                    current_match = self[i]
                    self.sort(False)
                    i = self.list.index(current_match)
                #i = self.list.index(overlap_matches[0]) - 1
            i = i+1
        self.sort(True)
        i=0
        finished = False
        while not finished:
            overlap = 0
            try:
                overlap = self[i].rend-self[i+1].rstart
                not_null = self[i].rstart!=self[i].rend and self[i+1].rstart!=self[i+1].rend
                not_the_same = self[i].rstart!=self[i+1].rstart and self[i].rend!=self[i+1].rend #ignore multicopy
                containment = self[i+1].rend<self[i].rend
            except IndexError:
                finished = True
                not_null = False
                containment = False
            if not_null and overlap>overlap_threshold:
                overlap_matches = self.find_opposite_overlaps(i, True)
                contain_overlap_1 = self.contain_interval(overlap_matches[0].qstart, overlap_matches[0].qend, False)
                for match in contain_overlap_1:
                    self.split_match(match, overlap_matches[0].qstart, overlap_matches[0].qend, False)
                contain_overlap_2 = self.contain_interval(overlap_matches[1].qstart, overlap_matches[1].qend, False)
                for match in contain_overlap_2:
                    self.split_match(match, overlap_matches[1].qstart, overlap_matches[1].qend, False)
                if containment:
                    current_match = self[i]
                    self.sort(True)
                    i = self.list.index(current_match)
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
    show_snps_output = subprocess.check_output(f"show-snps -TrH {prefix}.1delta", shell=True).strip().split(b'\n')

    assert(len(show_coords_output)>0)

    indels = []
    len_ref = -1
    len_query = -1
    if len(show_snps_output)>1:
        for line in show_snps_output:
            split_line = line.split(b'\t')
            rsub = split_line[1]
            qsub = split_line[2]
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
        og_matches.append(Match(start_ref, end_ref, start_query, end_query, strand_query, match_indels))

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
        jaccard_similarity = coverage_query/len_query
    else:
        jaccard_similarity = coverage_ref/len_ref
    jaccard_distance = 1-jaccard_similarity

    return plasmid_1_unimogs, plasmid_2_unimogs, jaccard_distance, blocks_ref, blocks_query

testing = True
if testing == True:
    #plasmid1 = "NZ_LT985234.1"
    #plasmid2 = "NZ_CP062902.1"
    #plasmid2 = "NZ_LR999867.1"
    #plasmid1 = "NZ_CP032890.1"
    #path = "/home/daria/Documents/projects/INC-plasmids/samples/fastas/incy"
    #plasmid2 = "cpe105_contig_5_np1212"
    #plasmid1 = "cpe061_contig_3_np1212"
    plasmid1 = "2_dup_3_dup"
    plasmid2 = "3_dup"
    path = "/home/daria/Documents/projects/pling/tests/test1/fastas"
    #path = "/home/daria/Documents/projects/mobmessing/plasmids_leah"
    print(new_integerise_plasmids(f"{path}/{plasmid1}.fna", f"{path}/{plasmid2}.fna", f"{plasmid1}~{plasmid2}", plasmid1, plasmid2, length_threshold=200))
    '''matches = Matches([Match(100578,102034,64267,65708,1),Match(101858,101881,57736,57759,1),Match(101881,102188,57759,58867,1)])
    matches.resolve_overlaps(0)
    length_threshold=15
    ref_to_block, query_to_block, max_id = make_interval_tree_w_dups(matches.list, length_threshold)
    populate_interval_tree_with_unmatched_blocks(ref_to_block, 102034, max_id+1, length_threshold)
    populate_interval_tree_with_unmatched_blocks(query_to_block, 65708, len(ref_to_block)+1, length_threshold)
    plasmid_1_unimogs = get_unimog(ref_to_block)
    plasmid_2_unimogs = get_unimog(query_to_block)
    print(plasmid_1_unimogs, plasmid_2_unimogs)'''
