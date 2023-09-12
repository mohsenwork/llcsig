import math
import numpy as np
import pytest
import networkx as nx
from simulations.data_generation.directed_graph import random_DAG, random_graph, strongly_connected_components, is_cyclic, extract_confs, random_edges, latent_edges, DirectedGraph


def test_latent_edges(nodes):
    n_obs, n_conf = nodes
    pairs = latent_edges(n_obs, n_conf)
    indicator = [0] * n_conf

    # Dimensions
    assert len(pairs) == 2 * n_conf

    for pair in pairs:

        indicator[pair[0]] += 1

        # test edge from latent to observed
        assert pair[0] in range(n_conf)
        assert pair[1] in range(n_conf, n_conf + n_obs)

    # Confounder have two children
    np.testing.assert_array_equal(indicator, [2] * n_conf)


def test_latent_edges_exception():
    n_obs = 100
    n_conf = math.comb(n_obs, 2) + 1

    # Too many confounders
    with pytest.raises(ValueError, match='Exceeded maximum number of latent pairs'):
        latent_edges(n_obs, n_conf)


def test_random_edges(nodes, lB, uB, deg):
    n_obs, n_conf = nodes
    node_list = list(range(n_conf, n_obs + n_conf))
    edges = random_edges(node_list, lB, uB, deg)
    degrees = {n: 0 for n in node_list}

    # no self cycles
    for edge in edges:

        # no self cycles
        assert edge[0] != edge[1]

        # edges between observed
        assert edge[0] in node_list
        assert edge[1] in node_list

        if deg:
            degrees[edge[1]] += 1
        else:
            degrees[edge[0]] += 1

    # correct degrees
    min_degree = min(degrees.values())
    max_degree = max(degrees.values())
    assert lB <= min_degree <= uB
    assert lB <= max_degree <= uB


def test_random_DAG(nodes, lB, uB):
    n_obs, n_conf = nodes
    n = n_obs + n_conf
    G = random_DAG(n_obs, n_conf, lB, uB)
    adj = nx.to_numpy_array(G).T

    # Acyclic
    assert not is_cyclic(G)
    # Dimensions
    np.testing.assert_array_equal(adj.shape, (n, n))
    # Bounds
    # Lower bound sometimes violated because graph is acyclic.
    assert np.max(np.count_nonzero(adj[n_conf:, n_conf:], axis=1)) <= uB


def test_random_graph(nodes, lB, uB, deg):
    n_obs, n_conf = nodes
    n = n_obs + n_conf
    G = random_graph(n_obs, n_conf, lB, uB, deg)
    adj = nx.to_numpy_array(G).T

    # No Self cycles
    np.testing.assert_array_equal(np.diag(adj), np.zeros(n))
    # Cyclic
    assert is_cyclic(G)
    # Dimensions
    np.testing.assert_array_equal(adj.shape, (n, n))
    # Bounds
    max_degree = np.max(np.count_nonzero(adj[n_conf:, n_conf:], axis=int(deg)))
    min_degree = np.min(np.count_nonzero(adj[n_conf:, n_conf:], axis=int(deg)))
    assert lB <= min_degree <= uB
    assert lB <= max_degree <= uB


def test_extract_confs(nodes):
    n_obs, n_conf = nodes
    edges_latent = latent_edges(n_obs, n_conf)
    # Graph with only confounders
    G = nx.DiGraph()
    G.add_nodes_from(range(n_obs + n_conf))
    G.add_edges_from(edges_latent)
    conf = extract_confs(G, n_obs, n_conf)

    # Dimensions
    np.testing.assert_array_equal(conf.shape, (n_obs, n_obs))

    # Content
    origins = {a for a, _ in edges_latent}
    for a in origins:
        dest = [j for i, j in edges_latent if i == a]

        n = dest[0] - n_conf
        m = dest[1] - n_conf
        assert conf[n, m] == 1
        assert conf[m, n] == 1


def test_init(nodes, lB, uB, cyclic, deg):
    n_obs, n_conf = nodes
    graph = DirectedGraph(n_obs, n_conf, cyclic, lB, uB, deg)

    # properties
    assert graph.n_conf == n_conf
    assert graph.n_obs == n_obs
    np.testing.assert_array_equal(graph.confs, extract_confs(graph._nx_graph, n_obs, n_conf))
    assert graph.cyclic == cyclic
    assert graph.cyclic == is_cyclic(graph._nx_graph)
    np.testing.assert_array_equal(graph.adjacency_observed, graph.adjacency[graph.n_conf:, graph.n_conf:, ])
    np.testing.assert_array_equal(graph.adjacency, nx.to_numpy_array(graph._nx_graph).T)


def test_strongly_connected_components_empty():
    # weak test since code depends on well tested package
    n = 10
    empty = np.zeros((n, n))
    true = [set([i]) for i in range(n)]
    pred = strongly_connected_components(empty)
    assert set(map(tuple, true)) == set(map(tuple, pred))


def test_strongly_connected_components_tril():
    # weak test since code depends on well tested package
    n = 10
    adj = np.tril(np.ones((10, 10)), k=-1)
    true = [set([i]) for i in range(n)]
    pred = strongly_connected_components(adj)
    assert set(map(tuple, true)) == set(map(tuple, pred))


def test_strongly_connected_components_case():
    # weak test since code depends on well tested package
    #   <-
    # 0 -> 1 -> 2
    adj = np.array([[0, 1, 0], [1, 0, 0], [0, 1, 0]])
    true = [{0, 1}, {2}]
    pred = strongly_connected_components(adj)

    assert set(map(tuple, true)) == set(map(tuple, pred))
