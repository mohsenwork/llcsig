# pylint: skip-file
from math import floor
import numpy as np
import pytest
from simulations.data_generation.linear_model import LinearModel, random_coefficients, random_coefficients_stable, is_stable
from simulations.data_generation.directed_graph import DirectedGraph, is_cyclic


def test_get_random_interventions(model: LinearModel):
    # Random interventions on 0 to n variables
    do = [model.get_random_interventions(i) for i in range(model.graph.n_obs + 1)]

    # Number of 1s
    np.testing.assert_array_equal(np.count_nonzero(do, axis=1), np.arange(model.graph.n_obs + 1))
    # Dimensions
    assert len(do[0]) == model.graph.n_obs
    assert len(do) == model.graph.n_obs + 1


def test_get_complete_interventions(model: LinearModel):
    do = model.get_complete_interventions()

    # Dimensions
    np.testing.assert_array_equal(do.shape, (model.graph.n_obs + 1, model.graph.n_obs))
    # Passive observatinal state
    assert 0 in do.sum(axis=1)
    # All single interventions
    np.testing.assert_array_equal(do.sum(axis=0), np.ones(model.graph.n_obs))


def test_get_indicator_observed(model: LinearModel):
    # Random interventions on half of variables
    n = model.graph.n_obs + model.graph.n_conf
    do = model.get_random_interventions(floor(model.graph.n_obs / 2))
    U = model._get_indicator_observed(do)
    diag = np.hstack((np.repeat(1, model.graph.n_conf), 1 - do))

    # Dimensions
    np.testing.assert_array_equal(U.shape, (n, n))
    # Content
    np.testing.assert_array_equal(diag, U.diagonal())


def test_get_error_matrix(model: LinearModel):
    n_samples = 1_000_000
    tol = 1  # one decimal
    n_intervened = floor(model.graph.n_obs / 2)

    # Random interventions on half of variables
    do = model.get_random_interventions(n_intervened)
    is_intervened = np.array(np.concatenate([[0] * model.graph.n_conf, do]), dtype=bool)
    errors = model._get_error_matrix(do, n_samples)

    # Intervened variables
    mean_intervened_sampled = errors[is_intervened].mean(axis=1)
    mean_intervened_true = np.repeat(0, n_intervened)
    cov_intervened_sampled = np.cov(errors[is_intervened])
    cov_intervened_true = np.diag(np.repeat(1, n_intervened))

    # Observed variables
    mean_obs_sampled = errors[np.logical_not(is_intervened)].mean(axis=1)
    mean_obs_true = model.disturbance_mean[np.logical_not(is_intervened)]
    cov_obs_sampled = np.cov(errors[np.logical_not(is_intervened)])
    cov_obs_true = np.diag(model.disturbance_std[np.logical_not(is_intervened)]**2)

    # Dimensions
    np.testing.assert_array_equal(errors.shape, (model.graph.n_conf + model.graph.n_obs, n_samples))

    # Means
    np.testing.assert_array_almost_equal(mean_intervened_sampled, mean_intervened_true, tol)
    np.testing.assert_array_almost_equal(mean_obs_sampled, mean_obs_true, tol)

    # Covariance
    np.testing.assert_array_almost_equal(cov_intervened_sampled, cov_intervened_true, tol)
    np.testing.assert_array_almost_equal(cov_obs_sampled, cov_obs_true, tol)


@pytest.mark.parametrize("n_obs, n_conf, low, high", [(5, 0, 0.1, 1.1), (50, 10, 1.3, 5.2)])
def test_random_coefficients(n_obs, n_conf, low, high):
    graph = DirectedGraph(n_obs, n_conf, True)
    coeffs = random_coefficients(graph.adjacency, low, high)
    _max = coeffs.max()
    _min = coeffs.min()

    # max
    assert (_max == 0) or (low <= _max <= high) or (-high <= _max <= -low)
    # min
    assert (_min == 0) or (low <= _min <= high) or (-high <= _min <= -low)
    # no edge = no coeff
    np.testing.assert_array_equal(graph.adjacency == 0, coeffs == 0)


@pytest.mark.parametrize("n", [10, 100, 1_000])
def test_is_stable_acyclic(n):
    graph = DirectedGraph(n, floor(n / 2), False)

    # Acyclic = stable
    assert is_stable(graph.adjacency)


def test_random_coefficients_stable_model(model: LinearModel):
    assert is_stable(model.coefficient_matrix)


def test_random_coefficients_stable(nodes):
    n_obs, n_conf = nodes
    graph = DirectedGraph(n_obs, n_conf, True)
    coeffs = random_coefficients_stable(graph.adjacency)
    assert is_stable(coeffs)


def test_sample(model: LinearModel):
    n_samples = 10**4
    do = model.get_random_interventions(floor(model.graph.n_obs / 2))
    samples = model.sample(n_samples, do)

    # Dimensions
    np.testing.assert_array_equal(samples.shape, (model.graph.n_obs, n_samples))


def test_get_model(model: LinearModel):
    n = model.graph.n_obs + model.graph.n_conf

    # Properties
    np.testing.assert_array_equal(model.graph.adjacency.shape, (n, n))
    np.testing.assert_array_equal(model.coefficient_matrix.shape, (n, n))
    np.testing.assert_array_equal(model.disturbance_mean.shape, n)
    np.testing.assert_array_equal(model.disturbance_mean.shape, n)
    assert is_cyclic(model.graph._nx_graph) == model.graph.cyclic
