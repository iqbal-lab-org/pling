from argparse import Namespace
from unittest import TestCase
import run_pling
from tests.utils import *


class Test_Pling_end_to_end(TestCase):
    def read_file(self, path):
        """Helper method to read a file."""
        with open(path, 'rb') as f:
            return f.read()

    def assert_files_are_identical(self, path_to_first_file, path_to_second_file):
        first_file_content = self.read_file(path_to_first_file)
        second_file_content = self.read_file(path_to_second_file)
        self.assertEqual(first_file_content, second_file_content)

    def test_pling_end_to_end(self):
        args = Namespace(genomes_list='tests/integration_test/data/incy_list_4.txt',
                         output_dir='tests/integration_test/data/out',
                         integerisation='align',
                         bakta_db=None,
                         jaccard=0.4,
                         dcj=4,
                         dedup=None,
                         dedup_threshold=None,
                         identity=80,
                         min_indel_size=200,
                         bh_connectivity=10,
                         bh_neighbours_edge_density=0.2,
                         small_subcommunity_size_threshold=4,
                         cores='2',
                         storetmp=False,
                         forceall=True,
                         ilp_solver='GLPK',
                         timelimit=None,
                         resources=None,
                         profile=None)
        run_pling.pling(args)

        assert_files_are_identical("tests/integration_test/data/out/all_plasmids_matrix.dist",
                                   "tests/integration_test/data/all_plasmids_matrix.truth.dist")
