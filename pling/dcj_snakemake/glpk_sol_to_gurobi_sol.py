def glpk_sol_to_gurobi_sol(input_fh, output_fh):
    print("# Solution for model obj", file=output_fh)

    gather_solution = False
    solution = []
    for line in input_fh:
        line = line.strip()

        is_objective_line = line.startswith("Objective:")
        if is_objective_line:
            obj_value = int(line.split()[3])
            print(f"# Objective value = {obj_value}", file=output_fh)

        solution_started = line == "No. Column name       Activity     Lower bound   Upper bound"
        if solution_started:
            gather_solution = True

        solution_ended = line.startswith("Integer feasibility conditions:")
        if solution_ended:
            gather_solution = False

        if gather_solution:
            solution.extend(line.split())

    headerless_solution = solution[13:]
    for number, column_name, star, activity, lower_bound, upper_bound in zip(*(iter(headerless_solution[i::6]) for i in range(0, 6))):
        print(column_name, activity, file=output_fh)


if __name__ == "__main__":
    import sys
    sys.exit(glpk_sol_to_gurobi_sol(sys.stdin, sys.stdout))
