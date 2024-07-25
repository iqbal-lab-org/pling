import subprocess

def unimog_to_ilp(unimog, lp, genome1, genome2):
    try:
        print("Converting unimog to ILP...\n")
        subprocess.run(f"dingII generate {unimog} -mm --writeilp {lp} -p {genome1} {genome2}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def ilp_GLPK(lp, solution, snakefile_dir, timelimit):
    try:
        print("Solving ILP with GLPK...\n")
        subprocess.run(f"glpsol --lp {lp} -o {solution}.tmp {timelimit}", shell=True, check=True, capture_output=True)
        print("Completed.\n Converting GLPK output to Gurobi output...\n")
        subprocess.run(f"python {snakefile_dir}/glpk_sol_to_gurobi_sol.py <{solution}.tmp >{solution}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
        subprocess.run(f"rm {solution}.tmp", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def dcj_dist(unimog, solution, genome1, genome2):
    try:
        print("Parsing Gurobi output to get DCJ distances...\n")
        dcj_out = subprocess.run(f"dingII parsesol {unimog} --solgur {solution} -p {genome1} {genome2}", shell=True, check=True, capture_output=True, text=True).stdout
        dist = int(dcj_out.strip().split(" ")[2])
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    return dist
