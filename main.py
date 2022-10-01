from collections import deque
import gurobipy as gp
from gurobipy import *
import bisect
import subprocess


class Graph:

    E = []
    V = []
    back_edges = {}
    forward_edges = {}

    def __init__(self, vertices, edges):

        # makes copies
        self.V = list(vertices)
        self.E = list(edges)

        self.back_edges = {}
        self.forward_edges = {}

        for v in self.V:
            self.back_edges[v] = []
            self.forward_edges[v] = []

        for e in self.E:
            # starting vertex is omitted from tuple
            self.forward_edges[e[0]].append((e[1], e[2], e[3], e[4]))
            self.back_edges[e[1]].append((e[0], e[2], e[3], e[4]))

    def N_forward(self, v):
        return self.forward_edges[v]

    def N_back(self, v):
        return self.back_edges[v]

    def get_V(self):
        return self.V

    def get_E(self):
        return self.E

    def __repr__(self):
        return "V = " + str(self.V) + "\nE =\n" + '\n'.join(map(str, self.E))


def find_in_sorted_list(elem, sorted_list):

    i = bisect.bisect_left(sorted_list, elem)
    if i != len(sorted_list) and sorted_list[i] == elem:
        return i
    return -1


def in_sorted_list(elem, sorted_list):

    i = bisect.bisect_left(sorted_list, elem)
    return i != len(sorted_list) and sorted_list[i] == elem


# obtain induced subgraph reachable from v with distance d
# V is resulting graph is not necessarily topologically sorted
def reachable_subgraph(G, v, d):

    active_nodes = deque()
    # (vertex, BFS level)
    active_nodes.append((v, 0))
    visited = {}
    V_ = [v]
    E_ = []

    # uses BFS to determine reachable subgraph in d steps
    while active_nodes:
        (u, level) = active_nodes.pop()
        for e in G.N_forward(u):
            E_.append((u, e[0], e[1], e[2], e[3]))

            if e[0] not in visited and level + 1 < d:
                active_nodes.append((e[0], level + 1))
                visited[e[0]] = 1
                V_.append(e[0])

    return get_induced_subgraph(G, V_)


def create_alignment_graph(G, S):

    V_a = []
    E_a = []
    for i in range(len(S) + 1):
        # add vertices for that level
        for v in G.get_V():
            V_a.append(v + "." + str(i))

        if i < len(S):
            # add insertion edges
            for e in G.get_E():
                E_a.append((e[0] + "." + str(i), e[1] + "." + str(i), 1, e[3], e[4]))

            # add deletion edges
            for v in G.get_V():
                E_a.append((v + "." + str(i), v + "." + str(i+1), 1, '-', e[4]))

            # add matching / sub edge
            for e in G.get_E():
                E_a.append((e[0] + "." + str(i), e[1] + "." + str(i+1), 0 if S[i] == e[2] else 1, e[3], e[4]))

    # add edges to sink
    V_a.append('end')
    for v in G.get_V():
        E_a.append((v + "." + str(len(S)), 'end', 0, '-', '-'))

    return Graph(V_a, E_a)


def create_reversed_graph(G):

    E_R = []
    for e in G.get_E():
        E_R.append((e[1], e[0], e[2], e[3], e[4]))
    return Graph(G.get_V(), E_R)


# remove multi-edges, keeping ones with lowest weight
def remove_multiedges(G):

    E_sorted = sorted(G.get_E())

    E_a_no_dup = [E_sorted[0]]
    for i in range(1, len(E_sorted)):
        e_prev = E_sorted[i-1]
        e = E_sorted[i]
        if e_prev[0] == e[0] and e_prev[1] == e[1]:
            continue
        else:
            E_a_no_dup.append(e)
    return Graph(G.get_V(), E_a_no_dup)


def get_induced_subgraph(G, V_s):

    V_s_sorted = sorted(V_s)
    E_s = []
    for e in G.get_E():
        # used sorted vertex list search to avoid quadratic check
        if in_sorted_list(e[0], V_s_sorted) and in_sorted_list(e[1], V_s_sorted):
            E_s.append(e)

    return Graph(V_s_sorted, E_s)


# returns topologically sorted list of the vertices in G
def topological_sort(G):

    counts = {}
    for v in G.get_V():
        counts[v] = len(G.N_back(v))

    sources = deque()
    for v in counts:
        if counts[v] == 0:
            sources.append(v)

    V_sorted = []
    while sources:
        u = sources.pop()
        V_sorted.append(u)
        for e in G.N_forward(u):
            counts[e[0]] = counts[e[0]] - 1
            if counts[e[0]] == 0:
                sources.append(e[0])

    return V_sorted


# returns a list of vertices reachable from the start vertex V[0] in at most delta steps
def vertices_reachable_from_start(G, delta):

    # topologically sort
    V_sorted = topological_sort(G)

    V_reachable = []
    min_distance = {}
    min_distance[V_sorted[0]] = 0
    for u in V_sorted:
        if u in min_distance and min_distance[u] <= delta:
            V_reachable.append(u)
        for e in G.N_forward(u):
            if e[0] not in min_distance:
                min_distance[e[0]] = min_distance[u] + e[1]
            else:
                min_distance[e[0]] = min(min_distance[e[0]], min_distance[u] + e[1])
    return V_reachable


def prune_alignment_graph(G, delta):

    # remove multi-edges keeping the ones with the lowest weight
    G_no_dup = remove_multiedges(G)

    # keep only vertices reachable from starting vertex with path of weight at most delta
    # and reachable from end in G^R with path of weight at most delta
    V_a_reachable_from_front = vertices_reachable_from_start(G_no_dup, delta)
    V_a_reachable_from_end = vertices_reachable_from_start(create_reversed_graph(G_no_dup), delta)

    # take intersection of the two lists
    V_a_pruned = list(set(V_a_reachable_from_front).intersection( V_a_reachable_from_end))

    return get_induced_subgraph(G, V_a_pruned)


def create_sub_ILP(model, G, delta, index_offset, global_var):

    # add ILP index to edges
    E = G.get_E()
    E_with_index = [(e[0], e[1], e[2], e[3], i + index_offset) for i, e in enumerate(E)]
    G_with_index = Graph(G.get_V(), E_with_index)

    # sort added to ensure V[0] is start vertex
    V = topological_sort(G_with_index)

    # Add sub-ILP variables
    y = model.addVars([e[4] for e in E_with_index], vtype=GRB.BINARY)

    # add constraints for source vertex
    rhs = gp.LinExpr(0)
    for e in G_with_index.N_forward(V[0]):
        rhs += y[e[3]]
    model.addConstr(1 == rhs)

    # add constraints for 'internal' vertices
    for i in range(1, len(V)):
        if V[i] != 'end':
            lhs = gp.LinExpr(0)
            for e in G_with_index.N_back(V[i]):
                lhs += y[e[3]]

            rhs = gp.LinExpr(0)
            for e in G_with_index.N_forward(V[i]):
                rhs += y[e[3]]
            model.addConstr(lhs == rhs)

    # add constraints for sink vertex
    lhs = gp.LinExpr(0)
    for e in G_with_index.N_back('end'):
        lhs += y[e[3]]
    model.addConstr(lhs == 1)

    # add weight constraint
    lhs = gp.LinExpr(0)
    for e in E_with_index:
        lhs += e[2] * y[e[4]]
    model.addConstr(lhs <= delta)

    # add global binding constraints
    for e in E_with_index:
        if e[3] != '-':
            lhs = gp.LinExpr(0)
            lhs += y[e[4]] + global_var[int(e[3])]
            model.addConstr(lhs <= 1)

    model.update()
    return model


# add extra component to edge tuple to store index for ILP variable
def add_component_to_E(E):

    E_added_compoenent = []
    for e in E:
        E_added_compoenent.append((e[0], e[1], e[2], e[3], '-'))
    return E_added_compoenent


# checks whether sample has variant at position pos
def variant_at_pos(pos, id, sample):

    vcf_file_name = 'chr_' + id + '_' + sample + '.vcf.gz'

    # GT=$($bcftools view -H chr_${id}_${sample}.vcf.gz | grep ${v}| cut -f 10)
    p1 = subprocess.Popen(['bcftools', 'view', '-H', vcf_file_name], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep', pos], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(['cut', '-f', '10'], stdin=p2.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p2.stdout.close()
    GT = p3.communicate()

    # TODO: incorporate phase information (for now we are just extracting the first phase)
    if GT[0] == '0':
        return False
    else:
        return True


# takes coordinate pos,  haplotype sample, length alpha, reference file name, reference id, length of reference
# returns substring from haplotype or empty string
def extract_substring(pos, sample, alpha, ref, id, ref_len):

    if not variant_at_pos(pos, id, sample):
        return ''
    else:
        vcf_file = 'chr_' + id + '_' + sample + '.vcf.gz'
        ref_file = ref + '.fa'

        # $samtools faidx ${ref}.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 1  chr${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa
        end_coordinate = str(min(int(pos) + alpha, ref_len))
        p1 = subprocess.Popen(['samtools', 'faidx', ref_file, id + ':' + pos + '-' + end_coordinate], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['bcftools', 'consensus', '-s', sample, 'H', '1', vcf_file],
                              stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output, err = p2.communicate()
        return output


# Takes graph G, location of variants L, alpha, delta, list of samples
def create_global_ILP(G, L, alpha, delta, number_variants, samples, ref, id, ref_len):

    model = gp.Model()
    global_var = model.addVars(range(1, number_variants+1), vtype=GRB.BINARY)
    index_offset = num_variants

    # add sub-ILPs for each variant position
    for pos in L:
        for sample in samples:
            S = extract_substring(pos, sample, alpha, ref, id, ref_len)
            if S != '':
                G_ind = reachable_subgraph(G, pos, alpha + delta)
                G_a = create_alignment_graph(G_ind, S)
                G_a_pruned = prune_alignment_graph(G_a, delta)
                model = create_sub_ILP(model, G_a_pruned, delta, index_offset, global_var)
                index_offset += len(G_a_pruned.get_E())

    # add global objective
    obj = gp.LinExpr(0)
    for i in range(1, number_variants+1):
        obj += global_var[i]
    model.setObjective(obj, GRB.MAXIMIZE)
    model.update()

    return model


if __name__ == "__main__":

    # vertices
    V = ['123', '124', '125', '126', '200']

    # edges, format: (start vertex, end vertex, symbol, variant #)
    E = [('123', '124', 'A', '-'), ('124', '125', 'T', '-'), ('123', '124', 'T', '1'), \
         ('123', '200', 'C', '2'), ('200', '125', 'G', '2'), ('125', '126', 'T', '-'), \
         ('125', '126', 'A', '3'), ('124', '126', 'T', '4')]

    # variant positions
    L = ['123', '124', '125']

    samples = ['HG00096']
    ref = 'hs37d5'
    id = '22'
    ref_len = 3
    num_variants = 4
    delta = 2
    alpha = 1

    # Program also assumes that vcf.gz files have been created for each sample in working directory
    # file name: chr_<id>_<sample>.vcf.gz

    # and reference file in working directory
    # file name: ref.fa

    # Need extra component in the edges for later ILP use
    # format: (start vertex, end vertex, symbol, variant #, ILP variable index)
    E = add_component_to_E(E)
    G = Graph(V, E)
    model = create_global_ILP(G, L, alpha, delta, num_variants, samples, ref, id, ref_len)
    print(model.display())
