from argparse import Namespace
from unittest import TestCase
from pling import run_pling
from tests.utils import *
from click.testing import CliRunner

class Test_Pling_end_to_end(TestCase):
    def read_file(self, path):
        """Helper method to read a file."""
        with open(path, "rb") as f:
            return f.read()

    def assert_files_are_identical(self, path_to_first_file, path_to_second_file):
        first_file_content = self.read_file(path_to_first_file)
        second_file_content = self.read_file(path_to_second_file)
        self.assertEqual(first_file_content, second_file_content)

    def test_pling_align_end_to_end(self):
        runner = CliRunner()

        result = runner.invoke(
            run_pling.cli,
            [
                "cluster", "align",
                "tests/integration_test/data/incy_list_4.txt",
                "tests/integration_test/data/out_align",
                "--containment_distance", "0.6",
                "--dcj", "4",
                "--cores", "2",
                "--forceall",
                "--batch_size", "50"
            ],
        )

        assert result.exit_code == 0, result.output

        assert_files_are_identical(
            "tests/integration_test/data/out_align/all_plasmids_distances.tsv",
            "tests/integration_test/data/all_plasmids_distances.align.truth.tsv",
        )

if __name__ == '__main__':
    unittest.main()
