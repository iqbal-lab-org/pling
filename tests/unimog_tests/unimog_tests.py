from argparse import Namespace
from unittest import TestCase
from pling import run_pling
from tests.utils import *
import os

class Test_Pling(TestCase):

    def setUp(self):
        self.toy_tests = ["test_1", "test_2"]
        self.toy_dir = "tests/unimog_tests/test_cases/toy_tests"

    def test_toy_alignment(self):
        tests = self.toy_tests
        for test in tests:
            args = Namespace(genomes_list=f"{self.toy_dir}/input/{test}.txt",
                             output_dir=f"{self.toy_dir}/{test}/out/align",
                             integerisation="align",
                             unimog=None,
                             containment_distance=0.2,
                             dcj=4,
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
                             sourmash_threshold=None,
                             topology=None,
                             regions=False
                             )
            run_pling.pling(args)

            assert_end_to_end(test, self.toy_dir)

class Test_Indels(TestCase):

    def setUp(self):
        self.indel_tests = ["test_1","test_2","test_3","test_4","test_5","test_6","test_7","test_8","test_9", "test_10"]
        self.indel_dir = "tests/unimog_tests/test_cases/indel_tests"

    def test_indel_alignment(self):
        tests = self.indel_tests
        for test in tests:
            args = Namespace(genomes_list=f"{self.indel_dir}/input/{test}.txt",
                             output_dir=f"{self.indel_dir}/{test}/out/align",
                             integerisation="align",
                             unimog=None,
                             containment_distance=0.6,
                             dcj=4,
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
                             sourmash_threshold=None,
                             topology=None,
                             regions=False
                             )
            run_pling.pling(args)

            #containment communities
            assert_containment(test, self.indel_dir, "align")
            #unimogs
            batch_num = 1
            assert_align_unimogs(test, self.indel_dir, batch_num)

class Test_Palindrome(TestCase):

    def setUp(self):
        self.test = "palindrome"
        self.dir = "tests/unimog_tests/test_cases/git_issues"

    def test_palindrome(self):
        args = Namespace(genomes_list=f"{self.dir}/input/{self.test}.txt",
                         output_dir=f"{self.dir}/{self.test}/out/align",
                         integerisation="align",
                         unimog=None,
                         containment_distance=0.6,
                         dcj=4,
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
                         sourmash_threshold=None,
                         topology=None,
                         regions=False
                         )
        run_pling.pling(args)

        assert_end_to_end(self.test, self.dir)

class Test_Containment_1(TestCase):

    def setUp(self):
        self.test = "containment_1"
        self.dir = "tests/unimog_tests/test_cases/git_issues"

    def test_containment_1(self):
        args = Namespace(genomes_list=f"tests/unimog_tests/test_cases/toy_tests/input/test_2.txt",
                         output_dir=f"{self.dir}/{self.test}/out/align",
                         integerisation="align",
                         unimog=None,
                         containment_distance=1,
                         dcj=4,
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
                         sourmash_threshold=None,
                         topology=None,
                         regions=False
                         )
        run_pling.pling(args)

        assert_end_to_end(self.test, self.dir)

if __name__ == '__main__':
    unittest.main()
