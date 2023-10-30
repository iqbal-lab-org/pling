from unittest import TestCase
from pling.dcj_snakemake.glpk_sol_to_gurobi_sol import glpk_sol_to_gurobi_sol
from tests.utils import *


class Test_glpk_sol_to_gurobi_sol(TestCase):
    def test_glpk_sol_to_gurobi_sol(self):
        with open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk.sol.txt") as input_fh, \
             open("tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk_to_gurobi.sol.txt", "w") as output_fh:
                glpk_sol_to_gurobi_sol(input_fh, output_fh)
        assert_files_are_identical(
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/glpk_to_gurobi.sol.txt",
            "tests/dcj_snakemake/test_cases/glpk_sol_to_gurobi_sol/gurobi.sol")
