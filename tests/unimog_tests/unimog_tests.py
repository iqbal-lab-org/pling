from argparse import Namespace
from unittest import TestCase
from pling import run_pling
from tests.utils import *
import os

def get_num_communities(filepath):
    with open(filepath) as communities_fh:
        num = len(communities_fh.readlines())
    return num

def get_batch_num(filepath):
    with open(filepath) as f:
        next(f)
        num = int(f.readline().strip())
    return num

def assert_containment(test, test_dir, integerisation):
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/all_pairs_containment_distance.tsv",
                               f"{test_dir}/{test}/truth/{integerisation}/jaccard/all_pairs_containment_distance.tsv")
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/containment_communities/objects/communities.tsv",
                               f"{test_dir}/{test}/truth/{integerisation}/jaccard/jaccard_communities/objects/communities.tsv")
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/containment_communities/objects/communities.txt",
                               f"{test_dir}/{test}/truth/{integerisation}/jaccard/jaccard_communities/objects/communities.txt")

def assert_anno_unimogs(test, test_dir, community_num, dedup):
    for community in range(community_num):
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/{community}_anno.unimog",
                                   f"{test_dir}/{test}/truth/anno/unimogs/{community}_anno.unimog")
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/{community}_map.txt",
                                   f"{test_dir}/{test}/truth/anno/unimogs/{community}_map.txt")
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks/{community}_blocks.unimog",
                                   f"{test_dir}/{test}/truth/anno/unimogs/relabelled/blocks/{community}_blocks.unimog")
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks/{community}_map_blocks.txt",
                                   f"{test_dir}/{test}/truth/anno/unimogs/relabelled/blocks/{community}_map_blocks.txt")
        if dedup:
            assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks_dedup/{community}_blocks.unimog",
                                       f"{test_dir}/{test}/truth/anno/unimogs/relabelled/blocks_dedup/{community}_blocks.unimog")
            assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks_dedup/{community}_map_blocks.txt",
                                       f"{test_dir}/{test}/truth/anno/unimogs/relabelled/blocks_dedup/{community}_map_blocks.txt")
            assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks_dedup/{community}_blocks.unimog",
                                       f"{test_dir}/{test}/truth/anno/unimogs/relabelled/dedup/{community}_dedup.unimog")
            assert_files_are_identical(f"{test_dir}/{test}/out/anno/unimogs/relabelled/blocks/{community}_map_dedup.txt",
                                       f"{test_dir}/{test}/truth/anno/unimogs/relabelled/dedup/{community}_map_dedup.txt")

def assert_align_unimogs(test, test_dir, batch_num):
    for i in range(batch_num):
        assert_files_are_identical(f"{test_dir}/{test}/out/align/unimogs/batch_{i}_align.unimog",
                                   f"{test_dir}/{test}/truth/align/unimogs/batch_{i}_align.unimog")
        assert_files_are_identical(f"{test_dir}/{test}/out/align/unimogs/batch_{i}_map.txt",
                                   f"{test_dir}/{test}/truth/align/unimogs/batch_{i}_map.txt")

def assert_annotation(tests, test_dir, dedup_bool):
    for test in tests:
        args = Namespace(genomes_list=f"{test_dir}/input/{test}.txt",
                         output_dir=f"{test_dir}/{test}/out/anno",
                         integerisation="anno",
                         bakta_db="tests/bakta_db/db-light",
                         containment_distance=0.2,
                         dcj=4,
                         dedup=dedup_bool,
                         dedup_threshold=None,
                         identity=80,
                         min_indel_size=200,
                         bh_connectivity=10,
                         bh_neighbours_edge_density=0.2,
                         small_subcommunity_size_threshold=4,
                         plasmid_metadata=None,
                         cores='2',
                         storetmp=False,
                         forceall=True,
                         ilp_solver="GLPK",
                         timelimit=None,
                         resources=None,
                         profile=None,
                         batch_size=5,
                         sourmash=False,
                         sourmash_threshold=None)
        run_pling.pling(args)

        #containment communities
        assert_containment(test, test_dir, "anno")
        #unimogs
        #community_num = get_num_communities(f"tests/{test}/truth/anno/containment/containment_communities.txt")
        #assert_anno_unimogs(test, community_num, dedup_bool)
        #DCJ distance matrix
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/all_plasmids_distances.tsv",
                                   f"{test_dir}/{test}/truth/anno/all_plasmids_distances.tsv")
        #DCJ communities
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/dcj_thresh_4_graph/objects/typing.tsv",
                                   f"{test_dir}/{test}/truth/anno/dcj_thresh_4_graph/objects/typing.tsv")
        assert_files_are_identical(f"{test_dir}/{test}/out/anno/dcj_thresh_4_graph/objects/hub_plasmids.csv",
                                   f"{test_dir}/{test}/truth/anno/dcj_thresh_4_graph/objects/hub_plasmids.csv")

class Test_Pling(TestCase):

    def setUp(self):
        self.toy_tests = ["test_1"]
        self.toy_dir = "tests/unimog_tests/test_cases/toy_tests"
        self.indel_tests = ["test_1","test_2","test_3","test_4","test_5","test_6","test_7","test_8","test_9"]
        self.indel_dir = "tests/unimog_tests/test_cases/indel_tests"

    '''
    def test_toy_annotation_no_dedup(self):
        assert_annotation([self.toy_tests[0]], self.toy_dir, False)
    '''

    def test_toy_alignment(self):
        tests = self.toy_tests
        for test in tests:
            args = Namespace(genomes_list=f"{self.toy_dir}/input/{test}.txt",
                             output_dir=f"{self.toy_dir}/{test}/out/align",
                             integerisation="align",
                             bakta_db=None,
                             containment_distance=0.2,
                             dcj=4,
                             dedup=None,
                             dedup_threshold=None,
                             identity=80,
                             min_indel_size=200,
                             bh_connectivity=10,
                             bh_neighbours_edge_density=0.2,
                             small_subcommunity_size_threshold=4,
                             plasmid_metadata=None,
                             cores='1',
                             storetmp=False,
                             forceall=True,
                             ilp_solver="GLPK",
                             timelimit=None,
                             resources=None,
                             profile=None,
                             batch_size=5,
                             sourmash=False,
                             sourmash_threshold=None)
            run_pling.pling(args)

            #containment communities
            assert_containment(test, self.toy_dir, "align")
            #unimogs
            batch_num = get_batch_num(f"{self.toy_dir}/{test}/out/align/batches/batching_info.txt")
            assert_align_unimogs(test, self.toy_dir, batch_num)
            #DCJ distance matrix
            assert_files_are_identical(f"{self.toy_dir}/{test}/out/align/all_plasmids_distances.tsv",
                                       f"{self.toy_dir}/{test}/truth/align/all_plasmids_distances.tsv")
            #DCJ communities
            assert_files_are_identical(f"{self.toy_dir}/{test}/out/align/dcj_thresh_4_graph/objects/typing.tsv",
                                       f"{self.toy_dir}/{test}/truth/align/dcj_thresh_4_graph/objects/typing.tsv")
            assert_files_are_identical(f"{self.toy_dir}/{test}/out/align/dcj_thresh_4_graph/objects/hub_plasmids.csv",
                                       f"{self.toy_dir}/{test}/truth/align/dcj_thresh_4_graph/objects/hub_plasmids.csv")

    def test_indel_alignment(self):
        tests = self.indel_tests
        for test in tests:
            args = Namespace(genomes_list=f"{self.indel_dir}/input/{test}.txt",
                             output_dir=f"{self.indel_dir}/{test}/out/align",
                             integerisation="align",
                             bakta_db=None,
                             containment_distance=0.6,
                             dcj=4,
                             dedup=None,
                             dedup_threshold=None,
                             identity=80,
                             min_indel_size=200,
                             bh_connectivity=10,
                             bh_neighbours_edge_density=0.2,
                             small_subcommunity_size_threshold=4,
                             plasmid_metadata=None,
                             cores='1',
                             storetmp=False,
                             forceall=True,
                             ilp_solver="GLPK",
                             timelimit=None,
                             resources=None,
                             profile=None,
                             batch_size=1,
                             sourmash=False,
                             sourmash_threshold=None)
            run_pling.pling(args)

            #containment communities
            assert_containment(test, self.indel_dir, "align")
            #unimogs
            batch_num = 1
            assert_align_unimogs(test, self.indel_dir, batch_num)


if __name__ == '__main__':
    unittest.main()
