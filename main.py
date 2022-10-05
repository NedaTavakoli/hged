import gurobipy as gp
import networkx as nx
import time
import sys


# obtain induced subgraph reachable from v with distance d
def reachable_subgraph(G, v, d):
    E = nx.bfs_edges(G, v, depth_limit=d)
    N = set([n for e in E for n in e])
    return nx.induced_subgraph(G, N)


def create_alignment_graph(G, start_v, S):

    G_A = nx.MultiDiGraph()
    for i in range(len(S) + 1):
        # add vertices for that level

        if i < len(S):
            # add insertion edges
            for e in G.edges:
                # start, end, (weight, variant)
                G_A.add_edge(e[0] + "." + str(i), e[1] + "." + str(i), (1, e[2][1]))

            # add deletion edges
            for v in G.nodes:
                G_A.add_edge(v + "." + str(i), v + "." + str(i+1), (1, '-'))

            # add matching / sub edge
            for e in G.edges:
                G_A.add_edge(e[0] + "." + str(i), e[1] + "." + str(i+1), (0 if S[i] == e[2][0] else 1, e[2][1]))

    # add edges to sink
    for v in G.nodes:
        G_A.add_edge(v + "." + str(len(S)), 'end', (0, '-'))

    return G_A, start_v + ".0", 'end'


# remove multi-edges, keeping ones with lowest weight
def remove_multiedges(G):

    E_sorted = sorted(G.edges)

    E_a_no_dup = [E_sorted[0]]
    for i in range(1, len(E_sorted)):
        e_prev = E_sorted[i-1]
        e = E_sorted[i]
        if e_prev[0] == e[0] and e_prev[1] == e[1]:
            continue
        else:
            E_a_no_dup.append(e)

    G_pruned = nx.MultiDiGraph()
    for e in E_a_no_dup:
        G_pruned.add_edge(e[0], e[1], e[2])

    return G_pruned


def prune_alignment_graph(G, start_x, end, delta):

    # remove multi-edges keeping the ones with the lowest weight
    G_no_dup = remove_multiedges(G)

    # keep only vertices reachable from starting vertex with path of weight at most delta
    # and reachable from end in G^R with path of weight at most delta
    # remove edges with weight 0
    path_lengths_from_front = nx.shortest_path_length(G_no_dup, source=start_x, weight=lambda _, __, a: list(a.keys())[0][0])
    path_lengths_from_end = nx.shortest_path_length(G_no_dup.reverse(), source=end, weight=lambda _, __, a: list(a.keys())[0][0])

    reachable_from_front = [v for v in path_lengths_from_front if path_lengths_from_front[v] <= delta]
    reachable_from_end = [v for v in path_lengths_from_end if path_lengths_from_end[v] <= delta]

    # take intersection of the two lists
    V_a_pruned = list(set(reachable_from_front).intersection(reachable_from_end))

    return nx.induced_subgraph(G_no_dup, V_a_pruned)


def create_sub_ILP(model, G, start_v, end_v, delta, index_offset, global_var):

    # add ILP index to edges
    E = G.edges
    G_with_index = nx.MultiDiGraph()
    for i, e in enumerate(E):
        # start, end, weight, variant, ILP-index
        G_with_index.add_edge(e[0], e[1], weight=e[2][0], variant=e[2][1], index=i+index_offset)

    # Add sub-ILP variables
    E_with_index = G_with_index.edges(data=True)
    y = model.addVars([e[2]['index'] for e in E_with_index], vtype=gp.GRB.BINARY)

    # add constraints for source vertex
    rhs = gp.LinExpr(0)
    for e in G_with_index.out_edges(start_v, data=True):
        rhs += y[e[2]['index']]
    model.addConstr(1 == rhs)

    # add constraints for 'internal' vertices
    V = G_with_index.nodes
    for i, v in enumerate(V):
        if v != 'end':
            lhs = gp.LinExpr(0)
            for e in G_with_index.in_edges(v, data=True):
                lhs += y[e[2]['index']]

            rhs = gp.LinExpr(0)
            for e in G_with_index.out_edges(v, data=True):
                rhs += y[e[2]['index']]
            model.addConstr(lhs == rhs)

    # add constraints for sink vertex
    lhs = gp.LinExpr(0)
    for e in G_with_index.in_edges(end_v, data=True):
        lhs += y[e[2]['index']]
    model.addConstr(lhs == 1)

    # add weight constraint
    lhs = gp.LinExpr(0)
    for e in G_with_index.edges(data=True):
        lhs += e[2]['weight'] * y[e[2]['index']]
    model.addConstr(lhs <= delta)

    # add global binding constraints
    for e in E_with_index:
        if e[2]['variant'] != '-':
            lhs = gp.LinExpr(0)
            lhs += y[e[2]['index']] + global_var[int(e[2]['variant'])]
            model.addConstr(lhs <= 1)

    model.update()
    return model


# Takes graph G, location of variants L, alpha, delta, list of samples
def create_global_ILP(G, locations, substrings, number_variants, alpha, delta):

    start = time.time()
    model = gp.Model()

    # add global variables
    global_var = model.addVars(range(1, number_variants+1), vtype=gp.GRB.BINARY)
    index_offset = num_variants

    # add sub-ILPs for each variant position
    for i, pos in enumerate(locations):
        for S in substrings[i]:

            print('\tadding sub-ILP for position number ' + str(i+1) + ' out of ' + str(len(locations)))
            start2 = time.time()
            G_ind = reachable_subgraph(G, pos, alpha + delta)
            end2 = time.time()
            total_eachable_subgraph = end2 - start2
            print('\tTotal time for reachable subgraph: ', total_eachable_subgraph)

            start3 = time.time()
            G_a, start_v, end_v = create_alignment_graph(G_ind, pos, S)
            end3 = time.time()
            total_create_alignment_graph = end3 - start3
            print('\tTotal time for creating alignment graph: ', total_create_alignment_graph)

            start4 = time.time()
            G_a_pruned = prune_alignment_graph(G_a, start_v, end_v, delta)
            end4 = time.time()
            total_prune_alignment_graph = end4 - start4
            print('\tTotal time for pruning alignment graph ', total_prune_alignment_graph)

            start5 = time.time()
            model = create_sub_ILP(model, G_a_pruned, start_v, end_v, delta, index_offset, global_var)
            end5 = time.time()
            total_create_sub_ILP = end5 - start5
            print('\tTotal time for create_sub_ILP: ', total_create_sub_ILP)

            index_offset += len(G_a_pruned.edges())

    # add global objective
    obj = gp.LinExpr(0)
    for i in range(1, number_variants+1):
        obj += global_var[i]
    model.setObjective(obj, gp.GRB.MAXIMIZE)
    model.update()

    # end time
    end = time.time()
    total_time_ILP = end - start
    print("Total time:", total_time_ILP)

    return model


def load_data(edges_file_name, location_substring_file_name):

    # load graph edges
    edge_file = open(edges_file_name, 'r')
    E = []
    for e in edge_file:
        e = e.strip().split()
        E.append((e[0], e[1], e[2], e[3]))
    edge_file.close()

    # load variant locations and substrings
    location_file = open(location_substring_file_name, 'r')
    locations = []
    substrings = []
    for l in location_file:
        arr = l.strip().split()
        locations.append(arr[0])
        substrings.append(arr[1:])
    location_file.close()

    # extract number of variants
    num_variants = len(set([e[3] for e in E if e[3] != '-']))

    return E, locations, substrings, num_variants


def construct_graph(E):
    G = nx.MultiDiGraph()
    for e in E:
        # start, end, (symbol, variant #)
        G.add_edge(e[0], e[1], (e[2], e[3]))
    return G

if __name__ == "__main__":


    # # edges, format: (start vertex, end vertex, symbol, variant #)
    # E = [('123', '124', 'A', '-'), ('124', '125', 'T', '-'), ('123', '124', 'T', '1'), \
    #      ('123', '200', 'C', '2'), ('200', '125', 'G', '2'), ('125', '126', 'T', '-'), \
    #      ('125', '126', 'A', '3'), ('124', '126', 'T', '4')]
    # locations = ['123']
    # substrings = [['TT']]
    # num_variants = 4
    # delta = 1
    # alpha = 2

    # arguments: edge_file, variant_location_substring_file, alpha, delta
    edges_file_name = sys.argv[1]
    location_substring_file_name = sys.argv[2]
    alpha = int(sys.argv[3])
    delta = int(sys.argv[4])
    print('Loading data...')
    E, locations, substrings, num_variants = load_data(edges_file_name, location_substring_file_name)

    print('Constructing graph...')
    G = construct_graph(E)

    print('Constructing ILP...')
    model = create_global_ILP(G, locations, substrings, num_variants, alpha, delta)
    print('Finished ILP construction')
    print('number variables: ' + str(model.NumVars))
    print('number constraints: ' + str(model.NumConstrs))

    print('\nSolving ILP...')
    model.optimize()
    num_variants_removed = sum(model.X[:num_variants])
    print('variants removed: ' + str(model.X[:num_variants]))
    print('Number variants removed: ' + str(num_variants_removed))
    print('Number variants retains: ' + str(num_variants - num_variants_removed))
