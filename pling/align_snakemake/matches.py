from intervaltree import IntervalTree, Interval
from pathlib import Path
import subprocess

class MatchPointsError(Exception):
    pass

class Indel:
    def __init__(self, rstart, qstart, len, type):
        self.type = type
        self.len = len
        self.rstart = rstart
        self.qstart = qstart
        if type == "DEL":
            self.rend = rstart+len
            self.qend = qstart+1
        elif type == "INS":
            self.rend = rstart+1
            self.qend = qstart+len

    def __str__(self):
        return f"({self.rstart}, {self.rend}, {self.qstart}, {self.qend}, {self.type})"

    def increase_len(self, len_increase):
        self.len += len_increase
        if self.type == "DEL":
            self.rend = self.rstart+self.len
        elif self.type == "INS":
            self.qend = self.qstart+self.len

class Match:
    def __init__(self, rstart, rend, qstart, qend, strand):
        if rend<rstart or qend<qstart:
            raise MatchPointsError(f"({rstart}, {rend}, {qstart}, {qend}, {strand}) is not a valid match! (start point greater than end point)")
        self.rstart = rstart
        self.rend = rend
        self.qstart = qstart
        self.qend = qend
        self.strand = strand

    def __setitem__(self, match):
        self = match

    def __eq__(self, other):
        if type(other) is type(self):
            return self.rstart == other.rstart and self.rend == other.rend and self.qstart == other.qstart and self.qend == other.qend
        return False

    def __hash__(self):
        return hash((self.rstart, self.rend, self.qstart, self.qend, self.strand))

    def __str__(self):
        return f"({self.rstart}, {self.rend}, {self.qstart}, {self.qend}, {self.strand})"

    def indel_intersect(self, indel):
        return ((self.rstart<indel.rend<=self.rend) or (self.rstart<=indel.rstart<self.rend)) and ((self.qstart<indel.qend<=self.qend) or (self.qstart<=indel.qstart<self.qend))

    def indel_at_rstart(self, indel):
        return indel.rstart<=self.rstart<indel.rend<=self.rend and ((self.qstart<indel.qend<=self.qend) or (self.qstart<=indel.qstart<=self.qend))

    def indel_at_qstart(self, indel):
        return indel.qstart<=self.qstart<indel.qend<=self.qend and ((self.rstart<indel.rend<=self.rend) or (self.rstart<=indel.rstart<=self.rend))

    def indel_at_rend(self, indel):
        return self.rstart<=indel.rstart<self.rend<=indel.rend and ((self.qstart<indel.qend<=self.qend) or (self.qstart<=indel.qstart<=self.qend))

    def indel_at_qend(self, indel):
        return self.qstart<=indel.qstart<self.qend<=indel.qend and ((self.rstart<indel.rend<=self.rend) or (self.rstart<=indel.rstart<=self.rend))

    def indel_strictly_contained(self, indel):
        return (self.rstart<=indel.rstart<=indel.rend<=self.rend) and (self.qstart<=indel.qstart<=indel.qend<=self.qend)

    def remove_indels_from_ends(self, indels):
        edited_bool = False
        for indel in indels:
            if self.indel_at_rstart(indel):
                edited_bool = True
                rstart = indel.rend
                if self.strand == 1:
                    qstart = indel.qend
                    edited = Match(rstart, self.rend, qstart, self.qend, self.strand)
                else:
                    qend = indel.qstart
                    edited = Match(rstart, self.rend, self.qstart, qend, self.strand)
            elif self.indel_at_qstart(indel):
                edited_bool = True
                qstart = indel.qend
                if self.strand == 1:
                    rstart = indel.rend
                    edited = Match(rstart, self.rend, qstart, self.qend, self.strand)
                else:
                    rend = indel.rstart
                    edited = Match(self.rstart, rend, qstart, self.qend, self.strand)
            elif self.indel_at_rend(indel):
                edited_bool = True
                rend = indel.rstart
                if self.strand == 1:
                    qend = indel.qend
                    edited = Match(self.rstart, rend, self.qstart, qend, self.strand)
                else:
                    qstart = indel.qend
                    edited = Match(self.rstart, rend, qstart, self.qend, self.strand)
            elif self.indel_at_qend(indel):
                edited_bool = True
                qend = indel.qstart
                if self.strand == 1:
                    rend = indel.rend
                    edited = Match(self.rstart, rend, self.qstart, qend, self.strand)
                else:
                    rstart = indel.rend
                    edited = Match(rstart, self.rend, self.qstart, qend, self.strand)
        if not edited_bool:
            edited = self
        return edited

    def projection(self, coord, indels, ref_bool): #if ref_bool=True, project from reference to query, else query to reference
        if ref_bool:
            if coord == self.rstart and self.strand == 1:
                return self.qstart
            elif coord == self.rstart and self.strand == -1:
                return self.qend
            elif coord == self.rend and self.strand == 1:
                return self.qend
            elif coord == self.rend and self.strand == -1:
                return self.qstart
        else:
            if coord == self.qstart and self.strand == 1:
                return self.rstart
            elif coord == self.qstart and self.strand == -1:
                return self.rend
            elif coord == self.qend and self.strand == 1:
                return self.rend
            elif coord == self.qend and self.strand == -1:
                return self.rstart
        if ref_bool:
            dist = coord - self.rstart
            for indel in indels:
                if self.indel_strictly_contained(indel):
                    if self.rstart<=indel.rstart<=coord<indel.rend:
                        return indel.qstart
                    elif self.rstart<=indel.rstart<coord:
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
            for indel in indels:
                if self.indel_strictly_contained(indel):
                    if self.qstart<=indel.qstart<=coord<indel.qend:
                        return indel.rstart
                    elif self.qstart<=indel.qstart<coord:
                        if indel.type == "INS":
                            dist -= indel.len
                        elif indel.type == "DEL":
                            dist += indel.len
            if self.strand == 1:
                projected_coord = self.rstart + dist
            else:
                projected_coord = self.rend - dist
        return projected_coord

class Matches:
    def __init__(self, list_of_matches, indels): #list_of_matches is a list of Match objects
        self.list = list_of_matches
        self.indels = indels
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

    def remove_indels_from_ends(self):
        for i in range(len(self.list)):
            self[i] = self[i].remove_indels_from_ends(self.indels)

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
        projected_start = match.projection(start, self.indels, ref_bool)
        projected_end = match.projection(end, self.indels, ref_bool)
        if ref_bool:
            if match.strand == 1:
                lhs_split = Match(match.rstart, start, match.qstart, projected_start, 1)
                interval = Match(start, end, projected_start, projected_end, 1)
                rhs_split = Match(end, match.rend, projected_end, match.qend, 1)
            else:
                rhs_split = Match(match.rstart, start, projected_start, match.qend, -1 )
                interval = Match(start, end, projected_end, projected_start, -1)
                lhs_split = Match(end, match.rend, match.qstart, projected_end, -1)
            if not lhs_split.qstart<=lhs_split.qend<=interval.qstart<=interval.qend<=rhs_split.qstart<=rhs_split.qend:
                raise Exception("While splitting match order was broken.")
        else:
            if match.strand == 1:
                lhs_split = Match(match.rstart, projected_start, match.qstart, start,  1)
                interval = Match(projected_start, projected_end, start, end, 1)
                rhs_split = Match(projected_end, match.rend, end, match.qend, 1)
            else:
                rhs_split = Match(projected_start, match.rend, match.qstart, start, -1 )
                interval = Match(projected_end, projected_start, start, end, -1)
                lhs_split = Match(match.rstart, projected_end, end, match.qend, -1)
            if not lhs_split.rstart<=lhs_split.rend<=interval.rstart<=interval.rend<=rhs_split.rstart<=rhs_split.rend:
                raise Exception("While splitting match order was broken.")
        if lhs_split.rend != lhs_split.rstart and lhs_split.qend!=lhs_split.qstart: #don't add interval
            lhs_split = lhs_split.remove_indels_from_ends(self.indels)
            interval = interval.remove_indels_from_ends(self.indels)
            self[index] = lhs_split
            self.insert(index+1, interval)
            if rhs_split.rend != rhs_split.rstart and rhs_split.qend!=rhs_split.qstart: #don't add null interval
                rhs_split = rhs_split.remove_indels_from_ends(self.indels)
                self.insert(index+2, rhs_split)
        else:
            interval = interval.remove_indels_from_ends(self.indels)
            self[index] = interval
            if rhs_split.rend != rhs_split.rstart and rhs_split.qend!=rhs_split.qstart: #don't add null interval
                rhs_split = rhs_split.remove_indels_from_ends(self.indels)
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
            projected_rstart_2 = self[i].projection(rstart_2, self.indels, ref_bool)
            projected_rend_1 = self[i].projection(rend_1, self.indels, ref_bool)
            if self[i].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, projected_rstart_2, projected_rend_1, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, projected_rend_1, projected_rstart_2, -1))
            projected_rend_1 = self[i+1].projection(rend_1, self.indels, ref_bool)
            projected_rstart_2 = self[i+1].projection(rstart_2, self.indels, ref_bool)
            if self[i+1].strand == 1:
                overlaps.append(Match(rstart_2, rend_1, projected_rstart_2, projected_rend_1, 1))
            else:
                overlaps.append(Match(rstart_2, rend_1, projected_rend_1, projected_rstart_2, -1))
        else:
            if qend_1>qend_2:
                qend_1 = self[i+1].qend
            projected_qstart_2 = self[i].projection(qstart_2, self.indels, ref_bool)
            projected_qend_1 = self[i].projection(qend_1, self.indels, ref_bool)
            if self[i].strand == 1:
                overlaps.append(Match(projected_qstart_2, projected_qend_1, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(projected_qend_1, projected_qstart_2, qstart_2, qend_1, -1))
            projected_qend_1 = self[i+1].projection(qend_1, self.indels, ref_bool)
            projected_qstart_2 = self[i+1].projection(qstart_2, self.indels, ref_bool)
            if self[i+1].strand == 1:
                overlaps.append(Match(projected_qstart_2, projected_qend_1, qstart_2, qend_1, 1))
            else:
                overlaps.append(Match(projected_qend_1, projected_qstart_2, qstart_2, qend_1, -1))
        return overlaps

    def resolve_overlaps(self, overlap_threshold):
        length = len(self.list)
        self.remove_indels_from_ends()
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
                try:
                    overlap_matches = self.find_opposite_overlaps(i, False)
                    contain_overlap_1 = self.contain_interval(overlap_matches[0].rstart, overlap_matches[0].rend, True)
                    for match in contain_overlap_1:
                        self.split_match(match, overlap_matches[0].rstart, overlap_matches[0].rend, True)
                    contain_overlap_2 = self.contain_interval(overlap_matches[1].rstart, overlap_matches[1].rend, True)
                    for match in contain_overlap_2:
                        self.split_match(match, overlap_matches[1].rstart, overlap_matches[1].rend, True)
                except MatchPointsError:
                    return
                if containment:
                    current_match = self[i]
                    self.sort(False)
                    i = self.list.index(current_match)
            i = i+1
            if i>length**3:
                raise Exception(f"The resolve_overlaps function in matches.py has been stuck in the while loop for {length**3} iterations! There is likely a bug, please raise an issue.")
        self.sort(True)
        i=0
        finished = False
        while not finished:
            overlap = 0
            try:
                overlap = self[i].rend-self[i+1].rstart
                not_null = self[i].rstart!=self[i].rend and self[i+1].rstart!=self[i+1].rend
                containment = self[i+1].rend<self[i].rend
            except IndexError:
                finished = True
                not_null = False
                containment = False
            if not_null and overlap>overlap_threshold:
                try:
                    overlap_matches = self.find_opposite_overlaps(i, True)
                    contain_overlap_1 = self.contain_interval(overlap_matches[0].qstart, overlap_matches[0].qend, False)
                    for match in contain_overlap_1:
                        self.split_match(match, overlap_matches[0].qstart, overlap_matches[0].qend, False)
                    contain_overlap_2 = self.contain_interval(overlap_matches[1].qstart, overlap_matches[1].qend, False)
                    for match in contain_overlap_2:
                        self.split_match(match, overlap_matches[1].qstart, overlap_matches[1].qend, False)
                except:
                    return
                if containment:
                    current_match = self[i]
                    self.sort(True)
                    i = self.list.index(current_match)
            i = i+1
            if i>length**3:
                raise Exception(f"The resolve_overlaps function in matches.py has been stuck in the while loop for {length**3} iterations! There is likely a bug, please raise an issue.")

testing = False
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
