from argparse import Namespace
from unittest import TestCase
from pling import run_pling
from tests.utils import *
from click.testing import CliRunner
import os

class Test_Pling(TestCase):

    def setUp(self):
        self.toy_tests = ["test_1"]
        self.toy_dir = "tests/unimog_tests/test_cases/toy_tests"

    def test_toy_alignment(self):
        runner = CliRunner()

        for test in self.toy_tests:
            result = runner.invoke(
                run_pling.cli,
                [
                    "cluster", "align",
                    f"{self.toy_dir}/input/{test}.txt",
                    f"{self.toy_dir}/{test}/out/align",
                    "--containment_distance", "0.2",
                    "--dcj", "4",
                    "--forceall",
                    "--batch_size", "5"
                ],
            )

            assert result.exit_code == 0, f"{test} failed:\n{result.output}"

            assert_end_to_end(test, self.toy_dir)

    def test_toy_skip(self):
        runner = CliRunner()

        for test in self.toy_tests:
            result = runner.invoke(
                run_pling.cli,
                [
                    "cluster", "skip",
                    f"{self.toy_dir}/input/{test}.txt",
                    f"{self.toy_dir}/{test}/out/skip",
                    f"{self.toy_dir}/input/skip.unimog",
                    "--containment_distance", "0.2",
                    "--dcj", "4",
                    "--forceall",
                    "--batch_size", "5"
                ],
            )

            assert result.exit_code == 0, f"{test} failed:\n{result.output}"

            dir = self.toy_dir

            #containment communities
            assert_containment(test, dir, "skip")

            #DCJ distance matrix
            assert_distances_identical(f"{dir}/{test}/out/skip/all_plasmids_distances.tsv",
                                        f"{dir}/{test}/truth/skip/all_plasmids_distances.tsv")
            #DCJ communities
            assert_files_are_identical(f"{dir}/{test}/out/skip/dcj_thresh_4_graph/objects/typing.tsv",
                                        f"{dir}/{test}/truth/skip/dcj_thresh_4_graph/objects/typing.tsv")
            assert_files_are_identical(f"{dir}/{test}/out/skip/dcj_thresh_4_graph/objects/hub_plasmids.csv",
                                        f"{dir}/{test}/truth/skip/dcj_thresh_4_graph/objects/hub_plasmids.csv")

class Test_Indels(TestCase):

    def setUp(self):
        self.indel_tests = ["test_1","test_2","test_3","test_4","test_5","test_6","test_7","test_8","test_9", "test_10"]
        self.indel_dir = "tests/unimog_tests/test_cases/indel_tests"

    def test_indel_alignment(self):
        runner = CliRunner()

        for test in self.indel_tests:
            print(test)
            result = runner.invoke(
                run_pling.cli,
                [
                    "cluster", "align",
                    f"{self.indel_dir}/input/{test}.txt",
                    f"{self.indel_dir}/{test}/out/align",
                    "--containment_distance", "0.6",
                    "--dcj", "4",
                    "--forceall",
                    "--batch_size", "1",
                ],
            )

            assert result.exit_code == 0, f"{test} failed:\n{result.output}"

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
        runner = CliRunner()

        result = runner.invoke(
            run_pling.cli,
            [
                "cluster", "align",
                f"{self.dir}/input/{self.test}.txt",
                f"{self.dir}/{self.test}/out/align",
                "--containment_distance", "0.6",
                "--dcj", "4",
                "--forceall",
                "--batch_size", "5",
            ],
        )

        assert result.exit_code == 0, result.output

        assert_end_to_end(self.test, self.dir)

class Test_Containment_1(TestCase):

    def setUp(self):
        self.test = "containment_1"
        self.dir = "tests/unimog_tests/test_cases/git_issues"

    def test_containment_1(self):
        runner = CliRunner()

        result = runner.invoke(
            run_pling.cli,
            [
                "cluster", "align",
                "tests/unimog_tests/test_cases/toy_tests/input/test_2.txt",
                f"{self.dir}/{self.test}/out/align",
                "--containment_distance", "1",
                "--dcj", "4",
                "--forceall",
                "--batch_size", "5",
            ],
        )

        assert result.exit_code == 0, result.output

        assert_end_to_end(self.test, self.dir)

if __name__ == '__main__':
    unittest.main()
