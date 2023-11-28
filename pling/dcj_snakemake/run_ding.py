import subprocess
from pling.utils import read_in_batch_pairs

def get_jaccard_distances_for_batch(jaccard_tsv):
    jaccards = {}
    with open(jaccard_tsv, "r") as f:
        for line in f:
            plasmid_1, plasmid_2, jaccard = line.strip().split("\t")
            jaccard = float(jaccard)
            jaccards[[plasmid_1,plasmid_2]] = jaccard
    return jaccard

def get_unimog(outputpath, integerisation, batch, genome1, genome2):
    unimog = ""
    if integerisation == "align":
        unimog = f"{OUTPUTPATH}/unimogs/batch_{batch}/{genome1}~{genome2}_align.unimog"
    '''
    elif integerisation == "anno":
        community = plasmid_to_community[genome1]
        unimog = f"{OUTPUTPATH}/unimogs/relabelled/blocks/{community}_blocks.unimog"
    '''
    return unimog

def unimog_to_ilp(unimog, lp, genome1, genome2, log):
    try:
        print("Converting unimog to ILP...\n")
        subprocess.run(f"dingII generate {unimog} -mm --writeilp {lp} -p {genome1} {genome2} 2>{log}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def ilp_gurobi(lp, solution, timelimit, log):
    try:
        print("Solving ILP with Gurobi...\n")
        subprocess.run(f"gurobi_cl ResultFile={solution} Threads=1 LogFile={log} {timelimit} {lp} >/dev/null", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def ilp_GLPK(lp, solution, snakefile_dir, timelimit, log):
    try:
        print("Solving ILP with GLPK...\n")
        subprocess.run(f"glpsol --lp {lp} -o {solution}.tmp {timelimit} 2>&1 > {log}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    try:
        print("Converting GLPK output to Gurobi output...\n")
        subprocess.run(f"python {snakefile_dir}/glpk_sol_to_gurobi_sol.py <{solution}.tmp >{solution}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    try:
        subprocess.run(f"rm {solution}.tmp", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def dcj_dist(unimog, solution, out_unimog_relabeled, genome1, genome2, dist_file, log):
    try:
        print("Parsing Gurobi output to get DCJ distances...\n")
        subprocess.run(f"dingII parsesol {unimog} --solgur {solution} --matching {out_unimog_relabeled} -p {genome1} {genome2} > {dist_file} 2>{log}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def batchwise_ding(pairs, jaccard_distance, jaccards, ilp_solver, integerisation, outputpath, batch, timelimit, snakefile_dir):
    for pair in pairs:
        genome1 = pair[0]
        genome2 = pair[1]
        lp = f"{outputpath}/tmp_files/ding/ilp/{genome1}~{genome2}.lp"
        solution = f"{outputpath}/tmp_files/ding/solutions/{genome1}~{genome2}.sol"
        dist_file = f"{outputpath}/tmp_files/dists_pairwise/{genome1}~{genome2}.dist"
        out_unimog_relabeled = f"{outputpath}/unimogs/matched/{genome1}~{genome2}_matched.unimog"
        log = "logs/unimog_to_ilp/{genome1}~{genome2}.log"
        if jaccards[pair]<=jaccard_distance:
            unimog = get_unimog(integerisation, batch, genome1, genome2)
                unimog_to_ilp(unimog, lp, genome1, genome2, log)
                if ilp_solver == "gurobi":
                    log = "logs/ilp/gurobi/{genome1}~{genome2}.log"
                    ilp_gurobi(lp, solution, timelimit, log)
                elif ilp_solver == "GLPK":
                    log = "logs/ilp/GLPK/{genome1}~{genome2}.log"
                    ilp_GLPK(lp, solution, snakefile_dir, timelimit, log)
                else:
                    raise RuntimeError(f"Unknown ILP solver: {ilp_solver}")
                log = "logs/dcj_dist/{genome1}~{genome2}.log"
                dcj_dist(unimog, solution, out_unimog_relabeled, genome1, genome2, dist_file, log)
    with open(f"{outputpath}/tmp_files/ding/completion/batch_{batch}", "w+") as f:
        f.write("done!")


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, jaccard and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--jaccard_tsv", required=True)
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--integerisation", required=True)
    parser.add_argument("--ilp_solver", required=True)
    parser.add_argument("--timelimit", required=True)
    parser.add_argument("--snakefile_dir", required=True)

    # Parse the arguments
    args = parser.parse_args()

    pairs=read_in_batch_pairs(f"{args.outputpath}/tmp_files/batches/batch_{args.batch}.txt")
    jaccards=get_jaccard_distances_for_batch()

    batchwise_ding(pairs, args.jaccard_distance, args.ilp_solver, args.integerisation, args.outputpath, args.batch, args.timelimit, args.snakefile_dir)
