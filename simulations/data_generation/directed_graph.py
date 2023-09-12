import itertools
import random
from copy import deepcopy
from typing import List, Set, Tuple
import math
import numpy as np
import networkx as nx
from networkx.exception import NetworkXNoCycle


def random_DAG(n_obs, n_conf, lower_bound: int, upper_bound: int) -> nx.DiGraph:
    '''
    Generates a random acyclic DiGraph.


    Return :
    -------
    An nx.Digraph object representing a directed acyclic graph.
    The graph has (n_obs + n_conf) columns and rows.
    The first n_conf rows are reserved for latent nodes and the rest for the observed.
    Latent nodes will always have no parents.

    Parameters:
    ----------
    n_obs : int
        The number of observed variables.
    n_conf : int
        The number of latent variable confounders.
    lower_bound : number
        The lower limit for the number of ingoing edges.
    upper_bound : number
        The upper limit for the number of ingoing edges.
    '''
    # random topological order
    topological_order = random.sample(range(n_conf, n_obs + n_conf), n_obs)
    pres = [topological_order[:i] for i, _ in enumerate(topological_order)]
    # random degree
    bounds = [(min([len(pre), lower_bound]), min([len(pre), upper_bound])) for pre in pres]
    degrees = [random.randint(lB, uB) for lB, uB in bounds]
    # random parents
    parents = [random.sample(pre, degree) for pre, degree in zip(pres, degrees)]
    # create edges
    edges_observed = [(p, v) for v, sublist in zip(topological_order, parents) for p in sublist]
    edges_latent = latent_edges(n_obs, n_conf)

    # create garph
    G = nx.DiGraph()
    G.add_nodes_from(range(n_obs + n_conf))
    G.add_edges_from(edges_observed + edges_latent)

    return G


def random_graph(n_obs: int, n_conf: int, lower_bound: int, upper_bound: int, indegree: bool) -> nx.DiGraph:
    '''
    Generates a random cyclic direct graph.


    Return :
    -------
    An nx.Digraph object representing a directed acyclic graph.
    A random graph respecting the bounds is created until it contains at least one cycle.
    The graph has (n_obs + n_conf) columns and rows.
    The first n_conf rows are reserved for latent nodes and the rest for the observed.
    Latent nodes will always have no parents.


    Parameters:
    ----------
    n_obs : int
        The number of observed variables.
    n_conf : int
        The number of latent variable confounders.
    lower_bound : number
        The lower limit for the number of outgoing or ingoing edges, depending on indegree.
    upper_bound : number
        The upper limit for the number of outgoing or ingoing edges, depending on indegree.
    indegree : bool
        If True, then the lower and upper bound parameters specify the indegree, else outdegree.
    '''
    while True:
        edges_observed = random_edges(range(n_conf, n_obs + n_conf), lower_bound, upper_bound, indegree)
        # create garph
        G = nx.DiGraph()
        G.add_nodes_from(range(n_obs + n_conf))
        G.add_edges_from(edges_observed)

        if is_cyclic(G):
            break

    edges_latent = latent_edges(n_obs, n_conf)
    G.add_edges_from(edges_latent)

    return G


def random_edges(nodes: List[int], lower_bound: int, upper_bound: int, indegree: bool) -> List[Tuple[int, int]]:
    '''
    Creates a list of random edges between nodes.
    A edge from a to b is represented by a tuple (a, b).
    Directed cycles are not allowed.


    Parameters:
    ----------
    nodes : List[int]
        A list or iterable of nodes.
    lower_bound : int
        The lower limit for the number of outgoing or ingoing edges, depending on indegree.
    upper_bound : int
        The upper limit for the number of outgoing or ingoing edges, depending on indegree.
    indegree : bool
        If True, then the lower and upper bound parameters specify the indegree, else outdegree.

    Return:
    ------
        A list of tuples.
    '''
    # random degree
    degrees = [random.randint(lower_bound, upper_bound) for _ in nodes]
    # random parents (no self cycles)
    candidates = [[el for el in nodes if el != node] for node in nodes]
    parents = [random.sample(sublist, degree) for degree, sublist in zip(degrees, candidates)]
    # create edges
    edges = [(p, v) for v, sublist in zip(nodes, parents) for p in sublist]

    # reverse edges for outdegree
    if not indegree:
        edges = [(v, p) for p, v in edges]

    return edges


def is_cyclic(G: nx.DiGraph) -> bool:
    '''
    Tests if a DiGraph is cyclic.


    Return:
    ------
    True if the graph contains a cycle, else False.

    Parameters:
    ----------
    G : nx.DiGraph
        A directed graph.
    '''
    try:
        nx.find_cycle(G)
    except NetworkXNoCycle:
        return False
    else:
        return True


def latent_edges(n_obs: int, n_conf: int) -> List[Tuple[int, int]]:
    '''
    Creates a list of tuples for edges from latent to observed.
    A tuple (a, b) represents an edge from a to b.
    Each latent variable has exactly two observed children.


    Return:
    ------
        A list of tuples.

    Parameters:
    ----------
    n_obs : int
        The number of observed variables.
    n_conf : int
        The number of latent variable confounders.
    '''
    if math.comb(n_obs, 2) < n_conf:
        raise ValueError('Exceeded maximum number of latent pairs')

    # select pairs of confounded nodes
    all_confs = list(itertools.combinations(range(n_conf, n_obs + n_conf), 2))
    np.random.shuffle(all_confs)
    confs = all_confs[:n_conf]

    # create list of edge tuples
    edges = [[(i, v1), (i, v2)] for i, (v1, v2) in enumerate(confs)]
    edges_flat = [item for sublist in edges for item in sublist]

    return edges_flat


def extract_confs(G: nx.DiGraph, n_obs: int, n_conf: int) -> np.ndarray:
    '''
    Extracts confounders between observed variables.


    Return:
    ------
    A n_obs x n_obs matrix G where G[i,j] = G[j,i] = 1 iff i and j are confounded.

    Parameters:
    ----------
    G : nx.DiGraph
        A directed graph.
    n_obs : int
        The number of observed variables.
    n_conf : int
        The number of latent variable confounders.
    '''
    conf = np.zeros((n_obs, n_obs))

    for i in range(n_conf):
        out_edges = [e[1] for e in G.edges(i)]
        v_1 = out_edges[0] - n_conf
        v_2 = out_edges[1] - n_conf
        conf[v_1, v_2] = 1
        conf[v_2, v_1] = 1

    return conf


def strongly_connected_components(adj: np.ndarray) -> List[Set[int]]:
    '''
    Extracts nodes in strongly connected components between nodes of an adjacency matrix.
    The adjacency matrix A has an entry A[i, j] = 1 iff there is an edge from j to i.


    Return:
    ------
    A list of sets where each set contains the indexes of a strongly connected component.
    '''
    G = nx.from_numpy_array(np.transpose(adj), create_using=nx.DiGraph)
    comps = nx.strongly_connected_components(G)
    return list(comps)


class DirectedGraph:

    @property
    def n_obs(self) -> int:
        # Number of observed variables
        return self._n_obs

    @property
    def n_conf(self) -> int:
        # Number of latent variables
        return self._n_conf

    @property
    def adjacency_observed(self) -> np.ndarray:
        # Adjacency matix of observed nodes
        # First observed variable has index 0
        # An entry i, j corresponds to an edge from j to i.
        return deepcopy(self._adjacency_observed)

    @property
    def confs(self) -> np.ndarray:
        # Confounder matix of observed nodes
        # Matrix where A[i,j] = A[j,i] = 1 iff i and j are confounded.
        return deepcopy(self._confs)

    @property
    def adjacency(self) -> np.ndarray:
        # Adjacency matix of observed and latent nodes
        # First confounder has index 0
        # An entry i, j corresponds to an edge from j to i.
        return deepcopy(self._adjacency)

    @property
    def cyclic(self) -> bool:
        # True if graph contains at least one cycle
        return self._cyclic

    def __init__(self, n_obs: int, n_conf: int, cyclic: int, lower_bound: int = 0,
                 upper_bound: int = 3, indegree: bool = True):
        '''
        Generates a random directed graph.


        Graph construction:
        ------------------
        The underlying graph is a networkx.DiGraph object, which has implementations of various useful methods
        such as serach methods and methods for computing graph properties (conncected components, etc..).
        If the graph is acyclic, a random topological order is created then edges are drawn according to the bounds.
        If the graph is cyclic, random edges are drawn according to the bounds until there is atleast one cycle.
        The first n_conf rows to are reserved for latent nodes and the rest for the observed.
        Self cycles are not allowed.

        Parameters:
        ----------
        n_obs : number
            The number of observed variables.
        n_conf : number
            The number of latent variable confounders. Each latent variable has exactly observed two children.
        cyclic : bool
            True if underlying graph should be cyclic.
        lower_bound : number
            The lower limit for the number of outgoing or ingoing edges, depending on indegree.
        upper_bound : number
            The upper limit for the number of outgoing or ingoing edges, depending on indegree.
        indegree : bool
            If True, then the lower and upper bound parameters specify the indegree, else outdegree.
        '''
        if cyclic:
            nx_graph = random_graph(n_obs, n_conf, lower_bound, upper_bound, indegree)
        else:
            nx_graph = random_DAG(n_obs, n_conf, lower_bound, upper_bound)

        self._cyclic = cyclic
        self._n_obs = n_obs
        self._n_conf = n_conf
        self._nx_graph = nx_graph
        self._adjacency = nx.to_numpy_array(nx_graph).T
        self._confs = extract_confs(nx_graph, n_obs, n_conf)
        obs_nodes = range(n_conf, n_obs + n_conf)
        self._adjacency_observed = nx.to_numpy_array(nx_graph.subgraph(obs_nodes)).T
