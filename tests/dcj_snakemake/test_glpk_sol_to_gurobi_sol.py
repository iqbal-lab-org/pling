from unittest import TestCase
from dcj_snakemake.glpk_sol_to_gurobi_sol import glpk_sol_to_gurobi_sol


class Test_glpk_sol_to_gurobi_sol(TestCase):
    def read_file(self, path):
        """Helper method to read a file."""
        with open(path, 'rb') as f:
            return f.read()

    def assert_files_are_identical(self, path_to_first_file, path_to_second_file):
        first_file_content = self.read_file(path_to_first_file)
        second_file_content = self.read_file(path_to_second_file)
        self.assertEqual(first_file_content, second_file_content)

    def test_glpk_sol_to_gurobi_sol(self):
        with open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk.sol.txt") as input_fh, \
             open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk_to_gurobi.sol.txt", "w") as output_fh:
                glpk_sol_to_gurobi_sol(input_fh, output_fh)
        self.assert_files_are_identical(
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk_to_gurobi.sol.txt",
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/gurobi.sol")
