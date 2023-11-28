from unittest import TestCase
from pling.dcj_snakemake.glpk_sol_to_gurobi_sol import glpk_sol_to_gurobi_sol
from tests.utils import *


class Test_glpk_sol_to_gurobi_sol(TestCase):
    def test_glpk_sol_to_gurobi_sol_test_case_short(self):
        with open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_1/glpk.sol.txt") as input_fh, \
             open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_1/glpk_to_gurobi.sol.txt", "w") as output_fh:
                glpk_sol_to_gurobi_sol(input_fh, output_fh)
        assert_files_are_identical(
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_1/glpk_to_gurobi.sol.txt",
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_1/gurobi.sol")

    def test_glpk_sol_to_gurobi_sol_test_case_long(self):
        with open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_2/glpk.sol.txt") as input_fh, \
             open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_2/glpk_to_gurobi.sol.txt", "w") as output_fh:
                glpk_sol_to_gurobi_sol(input_fh, output_fh)
        assert_files_are_identical(
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_2/glpk_to_gurobi.sol.txt",
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/test_case_2/gurobi.sol")