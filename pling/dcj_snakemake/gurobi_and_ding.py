import argparse
from pling.utils import read_in_batch_pairs
from pathlib import Path
from dingII.dingII_generate import *
from dingII.ilp_util_adj import *
from dingII.dingII_util import *
import gurobipy as gp
from gurobipy import GRB
from shared_functions import *

#modified function from ding to not take parser Namespace object as input
def read_unimog(unimog, genome1, genome2):
    with open(unimog) as f:
        lines = f.readlines()
    allgnms = readGenomes(lines)
    gnms = [g for g in allgnms if g[0] in [genome1, genome2]]
    if len(gnms) > 2:
        raise Exception('Multiple genomes of that name found!')
    elif len(gnms) < 2:
        raise Exception('Some genome names could not be found in the provided file!')
    return gnms

#modified function from ding to not take parser Namespace object as input
def knock_off_create_matching_model(gnms):
    amults = get_multiplicities(gnms[0])
    bmults = get_multiplicities(gnms[1])
    model = {}
    boundf = mm_bound
    model2 = create_model(amults, bmults, boundf)
    for gene in model2:
        if not gene in model:
            model[gene] = model2[gene]
    return model

#ripped off from dingII_generate main()
def make_ilp_dict(unimog, genome1, genome2):
    genomes = read_unimog(unimog, genome1, genome2)
    model = knock_off_create_matching_model(genomes)
    rd, gi1, gi2, circs, exts = full_relational_diagram(genomes, logging)
    logging.debug('Graph:')
    for c in rd.nodes(data=True):
        logging.debug(c)
    for c in rd.edges(data=True):
        logging.debug(c)
    sibs = siblings(rd, gi1, gi2)
    ilp = {}
    edge_constraints(rd, ilp)
    sibs_constraints(sibs, ilp, model)
    vertex_constraints(rd, ilp)
    singleton_constraints(circs, ilp)
    return ilp

def make_gurobi_model(unimog, genome1, genome2):
    ilp = make_ilp_dict(unimog, genome1, genome2)
    constraints = {k:v for k, v in ilp.items() if not k in ["binary", "general", "obj"]}
    m = gp.Model("ding")
    #add binary variables
    binary = [m.addVar(vtype=GRB.BINARY, name=var) for var in ilp["binary"]]
    m.update()
    binary_names = {binary[i].getAttr(GRB.Attr.VarName): (i,"b") for i in range(len(binary))}
    #add integer variables
    integer = [m.addVar(vtype=GRB.INTEGER, name=var[0], lb=var[1], ub=var[2]) for var in ilp["general"]]
    m.update()
    integer_names = {integer[i].getAttr(GRB.Attr.VarName): (i,"g") for i in range(len(integer))}
    names_dict = {**binary_names, **integer_names}
    #add constraints
    for key in constraints.keys():
        for constraint in constraints[key]:
            constraint_name = f"{key}.{constraint[0]}"
            og_elements = constraint[1].split(" ")
            elements = ['==' if el == '=' else el for el in og_elements]
            constants = [el for el in elements if el.isdigit()]
            var_names = [el for el in elements if not el.isdigit() and el not in ["==", "+", "-", "<", ">", "<=", ">="]]
            expression = ""
            for el in elements[:-1]:
                if el not in constants and el not in var_names:
                    expression += el
                elif el in var_names:
                    index=names_dict[el][0]
                    if names_dict[el][1]=="b":
                        var = f"binary[{index}] "
                    else:
                        var = f"integer[{index}] "
                    expression += var
                else:
                    expression += f"{el} * "
            expression += f" {elements[-1]}"
            m.addConstr(eval(expression), name = constraint_name)
    #add objective function
    obj = gp.LinExpr()
    for const, var_name in ilp["obj"]:
        if names_dict[var_name][1] == "b":
            var = binary[names_dict[var_name][0]]
        else:
            var = integer[names_dict[var_name][0]]
        obj += const * var
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    return m

def compute_DCJ(unimog, entry1, entry2, timelimit, threads):
    try:
        m = make_gurobi_model(unimog, entry1, entry2)
        if timelimit != None:
            m.setParam('TimeLimit', timelimit)
        m.setParam('Threads', threads)
        m.update()
        m.optimize()
        dist = m.ObjVal
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
        raise e
    except AttributeError:
        print('Encountered an attribute error')
        raise e
    return int(round(dist))

def batchwise_ding(pairs, jaccard_distance, jaccards, integerisation, outputpath, batch, timelimit, threads, plasmid_to_community):
    dists = []
    for pair in pairs:
        genome1 = pair[0]
        genome2 = pair[1]
        entry1, entry2 = get_entries(integerisation, genome1, genome2)
        if jaccards[(genome1,genome2)]<=jaccard_distance:
            unimog = get_unimog(outputpath, integerisation, plasmid_to_community, batch, genome1, genome2)
            dist = compute_DCJ(unimog, entry1, entry2, timelimit, threads)
            dists.append(f"{genome1}\t{genome2}\t{dist}\n")
        else:
            dists.append(f"{genome1}\t{genome2}\n")
    with open(f"{outputpath}/tmp_files/dists_batchwise/batch_{batch}_dcj.tsv", "w") as f:
        for line in dists:
            f.write(line)

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, jaccard and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--jaccard_tsv", required=True)
    parser.add_argument("--jaccard_distance", required=True)
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--communitypath", required=True)
    parser.add_argument("--integerisation", required=True, type=str)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--timelimit")

    # Parse the arguments
    args = parser.parse_args()

    pairs=read_in_batch_pairs(f"{args.outputpath}/batches/batch_{args.batch}.txt")
    jaccards=get_jaccard_distances_for_batch(args.jaccard_tsv)

    if args.integerisation=="anno":
        plasmid_to_community = get_plasmid_to_community(args.communitypath)
    else:
        plasmid_to_community=None

    batchwise_ding(pairs, float(args.jaccard_distance), jaccards, args.integerisation, args.outputpath, args.batch, args.timelimit, args.threads, plasmid_to_community)

if __name__ == "__main__":
    main()
