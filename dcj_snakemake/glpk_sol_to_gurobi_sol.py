import sys

print("# Solution for model obj")

parsing_solution = False
for line in sys.stdin:
    line = line.strip()

    is_objective_line = line.startswith("Objective:")
    if is_objective_line:
        obj_value = int(line.split()[3])
        print(f"# Objective value = {obj_value}")

    solution_parsing_started = line == "No. Column name       Activity     Lower bound   Upper bound"
    if solution_parsing_started:
        parsing_solution = True

    if parsing_solution:
        line_split = line.split()
        try:
            int(line_split[0])
        except (ValueError, IndexError):
            continue
        print(f"{line_split[1]} {line_split[3]}")
