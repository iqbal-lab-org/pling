def glpk_sol_to_gurobi_sol(input_fh, output_fh):
    print("# Solution for model obj", file=output_fh)

    parsing_solution = False
    for line in input_fh:
        line = line.strip()

        is_objective_line = line.startswith("Objective:")
        if is_objective_line:
            obj_value = int(line.split()[3])
            print(f"# Objective value = {obj_value}", file=output_fh)

        solution_parsing_started = line == "No. Column name       Activity     Lower bound   Upper bound"
        if solution_parsing_started:
            parsing_solution = True

        if parsing_solution:
            line_split = line.split()
            try:
                int(line_split[0])
            except (ValueError, IndexError):
                continue
            print(f"{line_split[1]} {line_split[3]}", file=output_fh)


if __name__ == "__main__":
    import sys
    sys.exit(glpk_sol_to_gurobi_sol(sys.stdin, sys.stdout))
