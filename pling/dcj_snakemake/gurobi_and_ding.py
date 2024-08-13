from pling.utils import read_in_batch_pairs
from pathlib import Path
from dingII.dingII_generate import *
from dingII.ilp_util_adj import *
from dingII.dingII_util import *
import gurobipy as gp
from gurobipy import GRB

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
            m.setParam('TimeLimit', float(timelimit))
        m.setParam('Threads', int(threads))
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
