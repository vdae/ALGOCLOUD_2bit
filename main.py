import networkx as nx
from os import listdir, path
import argparse
import math
from xml.dom.minidom import parse
import xml.dom.minidom
import numpy as np

INFINITY = float('inf')

def is_exact_k_connected(G, k):
    return nx.edge_connectivity(G) == k


def is_critical_k_connected(G, k):
    if k == 1:
        return nx.is_tree(G)
    if not is_exact_k_connected(G, k):
        return False
    for edge in G.edges:
        S = G.copy()
        S.remove_edge(edge[0], edge[1])
        if not is_exact_k_connected(G, k - 1):
            return False
    return True


def reduce_graph(Gr):
    G = Gr.copy()
    forbidden = []
    node = first_node_with_degree_2(G, forbidden)
    while node is not None:
        neighbors = [n for n in G.neighbors(node)]
        if len(neighbors) >= 2 and not G.has_edge(neighbors[0], neighbors[len(neighbors) - 1]):
            G.remove_node(node)
            G.add_edge(neighbors[0], neighbors[len(neighbors) - 1])
            node = first_node_with_degree_2(G, forbidden)
        else:
            forbidden.append(node)
            node = first_node_with_degree_2(G, forbidden)
    return G


def reduce_graph_multi(Gr):
    G = reduce_graph(Gr.copy())
    node = first_node_with_degree_2(G, [])
    while node is not None:
        neighbors = [n for n in G.neighbors(node)]
        if len(neighbors) >= 2:
            G.remove_node(node)
            G.add_edge(neighbors[0], neighbors[len(neighbors) - 1])
            node = first_node_with_degree_2(G, [])
        if reduced_finished(G):
            break
    return G


def reduced_finished(G):
    for (node, degree) in G.degree():
        neighbors = [n for n in G.neighbors(node)]
        if degree == 2 and len(neighbors) >= 2:
            return False
    return True


def first_node_with_degree_2(G, forbidden):
    for (node, degree) in G.degree():
        if degree == 2 and node not in forbidden:
            return node


def remove_nodes_deg_1(G):
    for (node, degree) in G.copy().degree():
        if degree == 1:
            G.remove_node(node)
    if 1 in [degree for (_, degree) in G.degree]:
        remove_nodes_deg_1(G)


def remove_nodes_deg_0(G):
    for (node, degree) in G.copy().degree():
        if degree == 0:
            G.remove_node(node)


def find_critical_pairs_multi(G):
    if first_node_with_degree_2(G, []) is not None:
        return ["not empty"]
    if nx.number_of_edges(G) <= 2:
        return ["not empty"]
    crit_pairs = []
    C = G.copy()
    for edge in C.edges:
        G.remove_edge(edge[0], edge[1])
        for other in C.edges:
            if edge == other:
                continue
            G.remove_edge(other[0], other[1])
            if not nx.is_connected(G):
                crit_pairs.append((edge, other))
            G.add_edge(other[0], other[1])
        G.add_edge(edge[0], edge[1])
    return crit_pairs


def is_outerplanar(G):
    OUT = nx.Graph(G)
    nodes = list(OUT.nodes)
    OUT.add_node("A")
    for node in nodes:
        OUT.add_edge("A", node)
    return nx.is_planar(OUT)


def remove_bridges(G):
    bridges = nx.bridges(G)
    for bridge in bridges:
        G.remove_edge(bridge[0], bridge[1])
    return [G.subgraph(c).copy() for c in nx.connected_components(G)]


def is_complete_bipartite(G):
    if nx.is_connected(G) and nx.is_bipartite(G):
        X, Y = nx.bipartite.sets(G)
        for x in X:
            for y in Y:
                if not G.has_edge(x, y):
                    return False
        return True
    return False

def find_resiliency(G):
    if is_outerplanar(G):
        return INFINITY
    elif (G.number_of_nodes() <= 5
          or nx.is_k_edge_connected(nx.Graph(G), 3)
          or is_complete_bipartite(G)
          or not find_critical_pairs_multi(G)):
        return 2
    else:
        return 1

def find_resiliency_reduced(G):
    if is_outerplanar(G):
        return INFINITY
    elif (nx.is_k_edge_connected(nx.Graph(G), 3)
          or not find_critical_pairs_multi(G)):
        return 2
    else:
        return 1


# SNDLib Parser by Da Silva et al.
# from https://github.com/carlosnatalino/python-simple-anycast-wdm-simulator
def read_sndlib_topology(file):
    graph = nx.Graph()

    with open('sndlib/' + file) as file:
        tree = xml.dom.minidom.parse(file)
        document = tree.documentElement

        graph.graph["coordinatesType"] = document.getElementsByTagName("nodes")[0].getAttribute("coordinatesType")

        nodes = document.getElementsByTagName("node")
        for node in nodes:
            x = node.getElementsByTagName("x")[0]
            y = node.getElementsByTagName("y")[0]
            graph.add_node(node.getAttribute("id"), pos=((float(x.childNodes[0].data), float(y.childNodes[0].data))))
        links = document.getElementsByTagName("link")
        for idx, link in enumerate(links):
            source = link.getElementsByTagName("source")[0]
            target = link.getElementsByTagName("target")[0]

            if graph.graph["coordinatesType"] == "geographical":
                length = np.around(calculate_geographical_distance(graph.nodes[source.childNodes[0].data]["pos"], graph.nodes[target.childNodes[0].data]["pos"]), 3)
            else:
                latlong1 = graph.nodes[source.childNodes[0].data]["pos"]
                latlong2 = graph.nodes[target.childNodes[0].data]["pos"]
                length = np.around(math.sqrt((latlong1[0] - latlong2[0]) ** 2 + (latlong1[1] - latlong2[1]) ** 2), 3)

            weight = 1.0
            graph.add_edge(source.childNodes[0].data, target.childNodes[0].data,
                           id=link.getAttribute("id"), weight=weight, length=length, index=idx)
    graph.graph["node_indices"] = []
    for idx, node in enumerate(graph.nodes()):
        graph.graph["node_indices"].append(node)

    return graph

# SNDLib Parser by Da Silva et al.
# from https://github.com/carlosnatalino/python-simple-anycast-wdm-simulator
def calculate_geographical_distance(latlong1, latlong2):
    R = 6373.0

    lat1 = math.radians(latlong1[0])
    lon1 = math.radians(latlong1[1])
    lat2 = math.radians(latlong2[0])
    lon2 = math.radians(latlong2[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    length = R * c
    return length

parser = argparse.ArgumentParser()
parser.add_argument("-pg", "--print-graphs", help="Print graphs and reduced graphs with matplotlib",
                    action=argparse.BooleanOptionalAction)
args = parser.parse_args()
num_graphs = 0
num_outer = 0
num_3con = 0
count = 0
count_bit = 0
print("TOPOLOGY ZOO: ")
for file in listdir("topology_zoo"):
    num_graphs += 1
    R = nx.read_graphml(path.join('topology_zoo', file))
    if not find_critical_pairs_multi(R):
        num_3con += 1
    if is_outerplanar(R):
        num_outer += 1
    remove_nodes_deg_1(R)
    remove_nodes_deg_0(R)
    if nx.is_empty(R):
        count += 1
        count_bit += 1
        continue
    remove_bridges(R)
    bi_components = list((R.subgraph(c).copy() for c in nx.biconnected_components(R)))
    connectivity = min(list(map(find_resiliency, bi_components)))
    if connectivity >= 2:
        count += 1
        count_bit += 1
        continue
    reduced_multi = list(map(reduce_graph_multi, bi_components))
    connectivity = min(list(map(find_resiliency_reduced, reduced_multi)))
    if connectivity >= 2:
        count_bit += 1
print("Number of Graphs: ", num_graphs)
print("Outerplanar: ", num_outer)
print("3-connected: ", num_3con)
print("Disassembled: ", count)
print("Reduced: ", count_bit)
print()

num_graphs = 0
num_outer = 0
num_3con = 0
count = 0
count_bit = 0
print("ROCKETFUEL: ")
for file in listdir("rf2zoo"):
    num_graphs += 1
    R = nx.read_gml(path.join('rf2zoo', file))
    if not find_critical_pairs_multi(R):
        num_3con += 1
    if is_outerplanar(R):
        num_outer += 1
    remove_nodes_deg_1(R)
    remove_nodes_deg_0(R)
    if nx.is_empty(R):
        count += 1
        count_bit += 1
        continue
    remove_bridges(R)
    bi_components = list((R.subgraph(c).copy() for c in nx.biconnected_components(R)))
    connectivity = min(list(map(find_resiliency, bi_components)))
    if connectivity >= 2:
        count += 1
    reduced_multi = list(map(reduce_graph_multi, bi_components))
    connectivity = min(list(map(find_resiliency, reduced_multi)))
    if connectivity >= 2:
        count_bit += 1
print("Number of Graphs: ", num_graphs)
print("Outerplanar: ", num_outer)
print("3-connected: ", num_3con)
print("Disassembled: ", count)
print("Reduced: ", count_bit)
print()

num_graphs = 0
num_outer = 0
num_3con = 0
count = 0
count_bit = 0
print("SNDLIB: ")
for file in listdir("sndlib"):
    num_graphs += 1
    R = read_sndlib_topology(file)
    if not find_critical_pairs_multi(R):
        num_3con += 1
    if is_outerplanar(R):
        num_outer += 1
    remove_nodes_deg_1(R)
    remove_nodes_deg_0(R)
    if nx.is_empty(R):
        count += 1
        count_bit += 1
        continue
    remove_bridges(R)
    bi_components = list((R.subgraph(c).copy() for c in nx.biconnected_components(R)))
    connectivity = min(list(map(find_resiliency, bi_components)))
    if connectivity >= 2:
        count += 1
    reduced_multi = list(map(reduce_graph_multi, bi_components))
    connectivity = min(list(map(find_resiliency, reduced_multi)))
    if connectivity >= 2:
        count_bit += 1
print("Number of Graphs: ", num_graphs)
print("Outerplanar: ", num_outer)
print("3-connected: ", num_3con)
print("Disassembled: ", count)
print("Reduced: ", count_bit)