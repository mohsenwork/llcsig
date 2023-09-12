# pylint: skip-file
import numpy as np
from pathlib import Path
from simulations.app import Simulations
from simulations.data_generation.linear_model import LinearModel
import shutil
import pickle
#import numpy as np


def test_against_llc_finite(model_params_small):
    # This test generate data from a model with very small noise then runs llc.
    # The output compared to the true model.
    n_obs, n_conf, cyclic, lower_bound, upper_bound, indegree, std = model_params_small
    n_samples = 1_000_000
    decimal = 5

    # create model and sample model
    model = LinearModel(n_obs, n_conf, cyclic, lower_bound, upper_bound, indegree, std)
    J = [[], [0], [1], [2], [3], [4]]
    data = Simulations.create_data(model, n_samples, J)

    # save model to tmp and run llc
    path = Path(Path.cwd(), 'tmp')
    path.mkdir(parents=True, exist_ok=True)
    Simulations.save_model(model, path)
    Simulations.run_llc(data=data, out_path=path, mode=1)

    # read llc output
    with Path(path, 'llc_out.pkl').open('rb') as f:
            llc_out = pickle.load(f)

    # remove tmp dir
    shutil.rmtree(path)

    # target and actual values
    B_est = np.array(llc_out['B'])
    B_true = model.coefficient_matrix[n_conf:, n_conf:]
    Confs_true = np.zeros((n_obs, n_obs))
    Confs_est = np.array(llc_out['Ce'])

    # Testing Coeffs	
    np.testing.assert_array_almost_equal(B_true, B_est, decimal=decimal)

    # Testing Confs
    np.testing.assert_array_almost_equal(Confs_true, Confs_est, decimal=decimal)
    
def test_against_llc_infinite(model):
    # This test generate the true covariance matrix of the observed nodes, and rns llc with it.
    # The output compared to the true model. LLC should find Σe and B when the pair and covariance 
    # conditions are met.
    
    decimal = 5
    single_interventions = model.get_complete_interventions()
    data = {tuple(interventions): model._get_covariance_observed_infinite(interventions)
            for interventions in single_interventions}

    # save model to tmp and run llc
    path = Path(Path.cwd(), 'tmp')
    path.mkdir(parents=True, exist_ok=True)
    Simulations.save_model(model, path)
    Simulations.run_llc(data=data, out_path=path, mode=3)
    
    # read llc output
    with Path(path, 'llc_out.pkl').open('rb') as f:
            llc_out = pickle.load(f)

    # remove tmp dir
    shutil.rmtree(path)

    # target and actual values
    B_est = np.array(llc_out['B'])
    B_true = model.coefficient_matrix[model.graph.n_conf:, model.graph.n_conf:]
    Confs_true = model._get_disturbance_covariance_matrix()
    Confs_est = np.array(llc_out['Ce'])

    # Testing Coeffs	
    np.testing.assert_array_almost_equal(B_true, B_est, decimal=decimal)

    # Testing Confs
    np.testing.assert_array_almost_equal(Confs_true, Confs_est, decimal=decimal)
