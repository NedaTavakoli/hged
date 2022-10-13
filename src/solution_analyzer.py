import sys


def get_variant_positions(edges_file_name):

    # load graph edges
    edge_file = open(edges_file_name, 'r')
    E = []
    end_of_ref = -1
    for e in edge_file:
        e = e.strip().split()
        E.append((e[0], e[1], e[2], e[3]))

        if e[3] == '-' and int(e[1]) > end_of_ref:
            end_of_ref = int(e[1])

    edge_file.close()

    positions = {}
    for e in E:
        if int(e[0]) <= int(end_of_ref) and e[3] != '-':
            positions[e[3]] = e[0]

    return positions


def get_retained_positions(ILP_sol_file_name, positions):

    ILP_sol_file = open(ILP_sol_file_name, 'r')
    sol_str = ILP_sol_file.read()
    ILP_sol_file.close()

    sol_str = sol_str.strip()
    sol_str = sol_str.replace('[', '')
    sol_str = sol_str.replace(']', '')
    sol_str = sol_str.replace(',', '')
    sol = map(int,map(float, sol_str.split()))

    retained_pos = []
    for i, x in enumerate(sol):
        if x == 0:
            retained_pos.append(positions[str(i)])

    return retained_pos


if __name__ == "__main__":

    edges_file_name = sys.argv[1]
    ILP_sol_file_name = sys.argv[2]
    output_file_name = sys.argv[3]

    print('Getting variant positions from graph')
    positions = get_variant_positions(edges_file_name)

    print('Getting retained variant positions from ILP solution')
    retained_positions = get_retained_positions(ILP_sol_file_name, positions)

    output_file = open(output_file_name, 'w')
    output_file.write('\n'.join(retained_positions))
    output_file.close()
