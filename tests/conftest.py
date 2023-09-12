# pylint: skip-file
from itertools import product
import logging
from collections import namedtuple
import pytest
from simulations.data_generation.directed_graph import DirectedGraph
from simulations.data_generation.linear_model import LinearModel

# Parameters for testing simulations/directed_graph.py, each configuration 
# consists of 4 parameters: n_obs, n_conf, high, low (input for functions random_graph() and random_DAG())

arr_nodes = [(5, 0), (5, 2), (20, 0), (20, 10), (50, 0), (50, 20)]
arr_cyclic = [True, False]
arr_lower_bound = [0, 1]
arr_upper_bound = [2, 3]
arr_indegree = [True, False]

@pytest.fixture(params=arr_nodes)
def nodes(request):
    return request.param

@pytest.fixture(params=arr_cyclic)
def cyclic(request):
    return request.param

@pytest.fixture(params=arr_lower_bound)
def lB(request):
    return request.param

@pytest.fixture(params=arr_upper_bound)
def uB(request):
    return request.param

@pytest.fixture(params=arr_indegree)
def deg(request):
    return request.param

# # Parameters for testing simulations/linear_model.py, each configuration 
# # consists of 7 parameters: n_obs, n_conf, cyclic, lower_bound, upper_bound, indegree, std (see linear_model constructor)
##  consider the cross product of all the configurations below for testing
# arr_obs_conf = [(5, 0), (5, 2), (20, 10), (50, 20)]
# arr_cyclic = [True, False]
# arr_lower_bound = [0, 1]
# arr_upper_bound = [3]
# arr_indegree = [True, False]
# arr_std = [1.5]
# model_params = [(*first, *rest) for first, *rest in product(arr_obs_conf,  arr_cyclic, arr_lower_bound, 
#                                                             arr_upper_bound, arr_indegree, arr_std)]
# # smaller test, above is more comprehensive but takes some time
model_params = [(5, 2, True, 0, 3, True, 1.5), (5, 2, False, 1, 2, False, 2)]

# Fixture for test of linear model
@pytest.fixture(params=model_params)
def model(request):
    yield LinearModel(*request.param)


# # Parameters for testing the entire code aganist the output of llc in the case of 
# # low error [standard diviation (1e-4)] and a large number of samples [1e6]
# # again, consider the cross product of the configurations for testing
# arr_obs_conf = [(5, 0), (5, 2)]
# arr_std = [1e-4, 1e-5]
# model_params_small = [(*first, *rest) for first, *rest in product(arr_obs_conf,  arr_cyclic, arr_lower_bound, 
#                                                                        arr_upper_bound, arr_indegree, arr_std)] 
# # # smaller test, above is more comprehensive but takes some time
model_params_small = [(5, 2, True, 0, 3, True, 1e-4)]#,(5, 0, True, 1, 2, True, 1e-5)]


# Fixture for test against llc output
@pytest.fixture(params=model_params_small)
def model_params_small(request):
    yield request.param
