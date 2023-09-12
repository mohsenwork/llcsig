import logging
import pickle
import random
import os
from typing import Dict, List, Tuple
import subprocess
from pathlib import Path
from timeit import default_timer as timer
import numpy as np
import pandas as pd
from .data_generation.linear_model import LinearModel
from .common.python.utilities import d_sep_CIRS, load_globals, pcor_dataframe, index_to_binary
from .resources.sigmasep.main import transform_tests_file_to_asp_file, run_all_CIT, run
logger = logging.getLogger(__name__)


class Simulations:

    @staticmethod
    def run(n_obs: int, n_conf: int, cyclic: bool, lB: int, uB: int, n_samples: List[int], experiments: Dict, models: List[int],
            name: str, asp_modes: List[int], seps: List[str], llc_mode: int, faithful_modes: List[int], penalty: str = 'none',
            reg: int = 0, n_bootstraps: int = 30, maxsetsize: int = -1, rules: List[int] = [1,2,3], null: int = 1, 
            alpha_llc: float = 0.05, alpha_asp: float = 0.01):
        '''
        Entry point for simulating and running LLC / ASP algorithms. 
        Combines model and data generation, runs the specified algorithms and saves the result.
        For an example check out __main__.py.


        Parameters:
        ----------
        n_obs: int,
            The number of observed variables.
        n_conf: int,
            The number of pair confounders.
        cyclic: bool,
            Whether the model should be cyclic or not.
        lB: int, uB: int,
            Lower and upper bound for the number of edges per node. 
            For more details, see linear_model.py and directed_graph.py.
        n_samples: List[int],
            List of sample numbers. Use np.inf to indicate the infinte sample limit.
            Example: [1000, np.inf]
        experiments: Dict,
            A dictionary specifying the experiments/interventions to be performed.
            The key are experiments names (data will be saved in a folder under that name)
            The values are a list of list of indexes of variables to intervene on:
            Example {'foo': [[], [0, 1], [2]]} stands for an experimental setup where Data is generated from three interventions:
            []: Data is generated from the passive observational state (zero intervention)
            [0, 1]: Data is generated from an intervention on variables 0 and 1
            [2]: Data is generated from an intervention on variable 2
            The algorithms will run and the output will be saved in a folder names 'experiment_foo'
        models: List[int],
            A list of model indexes, specifying how many models to generate.
            A model with index i will be saved under a folder containing that index.
        name: str,
            Name of the folder under which all data will be saved.
        asp_modes: List[int],
            A list of modes to run the ASP approach in.
            1: Run clingo to find one optimal answer set.
            2: Run clingo once to compute the union all optimal answer sets and once to compute their intersection.
            3: Run clingo once per possible features with additional hard constraints.
        seps: List[str],
            A list specifying which criterion(s) to use.
            's' corresponds to sigma-separantion encoding in sigma_hej_cyclic.pl, 
            'd' corresponds to d-separantion encoding for cyclic case in hej_cyclic.pl,
            'd' corresponds to d-separantion encoding for acyclic case in hej_acyclic.pl,
        llc_mode: int,
            Mode to run LLC in.
            mode=1: plain version with finite samples. Runs LLC once. The argument mode should be 1
                    The data in data_path must is pickled python dictionary with the intervention vector as keys 
                    and sampls as values. The covariance matrix is calculated, then llc_plain.R is called.
            mode=2: bootstrapped version. Runs LLC with resampling. The argument mode should be 2.
                    The data in data_path must is pickled python dictionary with the intervention vector as keys 
                    and sampls as values. llc_bootstrap.R is called with the samples. 
            mode=3: plain version with oracle covariance. Runs LLC once. The argument mode should be 3
                    The data in data_path must is pickled python dictionary with the intervention vector as keys 
                    and oracle cov matrixes as values. llc_plain.R is called using the oracle covariance matrix.
        faithful_modes: List[int],
            A list specifying if faithfulness rules should be added.
            0: Run LLC with no additional faithfulness rules.
            1: Run LLC with additional faithfulness rules specified by the rules parameter.
        penalty: str,
             Type of penalization/regularization used in the solving of the linear equation systems of LLC.
             Can be 'none', 'L0', 'L1' or 'L2'.
        reg: int,
            The regularization parameter used in the solving of the linear equation systems of LLC.
        n_bootstraps : int
            Number times to run llc_plain (only relevant for mode 2).
        maxsetsize : int
            The maximum size of a set that is used to condition on when using faithfulness rules. -1 for no restrictions.
        rules : List[int]
            Which of the faithfulness rules should be used. Any sublist of [1,2,3] can be used. An empty list means no faithfulness assumptions.
        null : int
            Specifies whether to request identifiability information about the null-space. Use null=1 to request such information (only relevant for modes 1 and 3).
        alpha_llc : float
            The significance level for the conditional independence tests when running LLC, used for faithfulness rules.
        alpha_asp : float
            The significance level for the conditional independence tests for ASP.
        '''
        home_path = Path.cwd()
        study_path = Path(home_path, 'simulations', 'data', f'{name}')
        asp_path = Path(home_path, 'simulations', 'resources', 'sigmasep', 'ASP')

        # load utility R code
        load_globals(Path(home_path, 'simulations', 'common', 'R', 'utilities.R'))

        # create main study directory
        study_path.mkdir(exist_ok=True, parents=True)
        logger.info('Starting simulations')

        for i in models:
            time_rows = []
            # set seeds
            np.random.seed(i)
            random.seed(i)

            # create model and save to file
            logger.info('Creating model %d ...', i)
            model = LinearModel(n_obs=n_obs, n_conf=n_conf, cyclic=cyclic, lower_bound=lB, upper_bound=uB)
            Simulations.save_model(model, Path(study_path, f'model_{i:03d}'))

            # random seed for sampling model
            np.random.seed()
            random.seed()

            # sample model
            logger.info('Sampling model %d ...', i)
            params = [(n, experiments[exp_name], exp_name) for n in n_samples for exp_name in experiments.keys()]

            # like (1000, 'name', [[], [0], [1]])
            for n, exp_list, exp_name in params:

                if n == np.inf or n == 'inf':

                    if len(seps) > 0:
                        # ASP input is conditional independence relations from d-separation
                        input_asp = d_sep_CIRS(model, exp_list)
                        use_CIT = False

                    # LLC input is oracle covariance matrix of observed nodes
                    exp_list_binary = index_to_binary(exp_list, model.graph.n_obs)
                    input_llc = {tuple(exp): model._get_covariance_observed_infinite(exp) for exp in exp_list_binary}
                    # Extract and save true partial correlation coefficients from covariance matrix
                    df = pcor_dataframe(input_llc)
                    p = Simulations.create_dir(study_path, i, n, exp_name)
                    df.to_csv(Path(p, f'true_p_cor.csv'))

                    mode_llc=3

                else:
                    # Sample model
                    data = Simulations.create_data(model, n, exp_list) # format the same; (0,0,0,1): samples
                    # ASP input is conditional independence relations from data
                    Data = [(np.transpose(S).copy(), np.nonzero(do)[0]) for do, S in data.items()]
                    if len(seps) > 0:
                        input_asp = run_all_CIT(Data, n_obs)
                    use_CIT = True
                    # LLC input is samples data
                    input_llc = data
                    mode_llc = llc_mode
                    
                # run llc without faithfulness rules
                if 0 in faithful_modes:
                    p = Simulations.create_dir(study_path, i, n, exp_name, method_name=f'llc_nf_{mode_llc}')
                    t = Simulations.run_llc(data=input_llc, out_path=p, mode=mode_llc, penalty=penalty, reg=reg, n_bootstraps=n_bootstraps, null=null, rules=[], alpha=alpha_llc)
                    time_rows += [{'method': f'llc_nf_{mode_llc}', 'n': n, 't': t, 'i': i, 'e': exp_name}]

                # run llc with given faithfulness rules
                if 1 in faithful_modes:  
                    p = Simulations.create_dir(study_path, i, n, exp_name, method_name=f'llc_f_{mode_llc}')
                    t = Simulations.run_llc(data=input_llc, out_path=p, mode=mode_llc, penalty=penalty, reg=reg, n_bootstraps=n_bootstraps, null=null, maxsetsize=maxsetsize, rules=rules, alpha=alpha_llc)
                    time_rows += [{'method': f'llc_f_{mode_llc}', 'n': n, 't': t, 'i': i, 'e': exp_name}]

                # Run ASP
                if len(seps) > 0:
                    modes_paths = [(mode, Simulations.create_dir(study_path, i, n, exp_name, f'asp_{mode}')) for mode in asp_modes]
                    time_rows += [{'method': f'asp_{mode}_{dic["sep"]}', 'n': n, 't': dic['t'], 'i': i, 'e': exp_name}
                                for mode, path in modes_paths
                                for dic in Simulations.run_asp(input_asp, path, asp_path, seps, n_obs, use_CIT, mode=mode, alpha=alpha_asp)]                            
                    
                    # Save CIRs
                    p = Simulations.create_dir(study_path, i, n, exp_name)
                    input_asp.to_csv(Path(p, f'cir_{"sep" if n == np.inf else "cit"}.csv'))

            # save times
            df = pd.DataFrame(time_rows)
            output_path = Path(study_path, 'times.csv')
            df.to_csv(output_path, mode='a', header=not os.path.exists(output_path))
            logger.info('Done.')

    @staticmethod
    def create_dir(study_path: Path, model_nr: int, n: int, exp_name: str, method_name: str = None):
        '''
        Creates data directory.


        Return:
        ------
        A paths to the created directory.

        Parameters:
        ----------
        study_path : Path
            Where to create data directories.
        model_nr : int
            Model number.
        n : int
            Number of samples.
        exp_name : str
            Experiment folder name.
        method_name : str
            Method folder name.
        '''
        if method_name is not None:
            path = Path(study_path, f'model_{model_nr:03d}', f'samplesize_{n}', f'experiment_{exp_name}', method_name)
        else:
            path = Path(study_path, f'model_{model_nr:03d}', f'samplesize_{n}', f'experiment_{exp_name}')
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def run_llc(data, out_path: Path, mode: int = 1, penalty: str = 'none', reg: int = 0, n_bootstraps: int = 30,
                maxsetsize: int = -1, rules: List[int] = [], null: int = 1, alpha: float = 0.05, verbose: int = 0) -> float:
        '''
        Runs the LLC algorithm.


        Return:
        ------
        LLC runtime.

        Parameters:
        ----------
        data : Dict[Tuple[int, ...] np.ndarray]
            Input data. The keys are tuples indicating the performed interventions and the values are either samples or a covariance matrix (see mode).
        out_path : Path
            Path to directory were to save the output.
        mode : int
            1 means run plain version of LLC with finite data. This runs LLC once. The values of the data dictionary are samples.
            2 means run bootstrapped version of LLC. This runs LLC multiple times with resampling. The values of the data dictionary are samples.
            3 means run plain version of LLC in infinite sample limit. This runs LLC once with oracle covariance matrices (values of the data dictionary).
        penalty : str
            'none', 'L0', 'L1' or 'L2' - type of penalization/regularization used in the solving of the linear equation systems.
        reg : int
            Regularization parameter.
        n_bootstraps : int
            Number times to run llc_plain (only relevant for mode 2).
        maxsetsize : int
            The maximum size of a set that is used to condition on when using faithfulness rules. -1 for no restrictions.
        rules : List[int]
            Which of the faithfulness rules should be used. Any sublist of [1,2,3] can be used. An empty list means no faithfulness assumptions.
        null : int
            Specifies whether to request identifiability information about the null-space. Use null=1 to request such information (only relevant for modes 1 and 3).
        alpha : float
            The significance level for the conditional independence tests, used for faithfulness rules.
        verbose : int
            Print updates to terminal iff verbose = 1.
        '''
        logger.info('Running LLC ...')
        data_path = Simulations.save_llc(data, out_path)
        start = timer()
        if len(rules) > 0:
            rules_str = ''.join([str(el) for el in rules])
        else:
            rules_str = '0'
        command = f'Rscript simulations/resources/llc/2012/main.R {data_path} {out_path.absolute()} {mode} {penalty} {reg} {n_bootstraps} {maxsetsize} {rules_str} {null} {alpha} {verbose}'
        subprocess.call(command, shell=True)
        t = timer() - start
        data_path.unlink()
        return t

    @staticmethod
    def run_asp(data, out_path: Path, asp_path: Path, seps: List[str], n_obs: int, use_CIT: bool, mode: int, alpha: float) -> List[float]:
        '''
        Runs the ASP algorithm.


        Return:
        ------
        Arrray of runtimes, one value for each element in seps.

        Parameters:
        ----------
        data : Dict[Tuple[int, ...] np.ndarray]
            Input data. The keys are tuples indicating the performed interventions and the values are the samples.
        out_path : Path
            Path to directory were to save the output.
        asp_path : Path
            Path were to ASP encoding directory.
        seps : List[str]
            List of encodings to query the ASP with.
        n_obs : int
            The number of observed variables.
        model : LinearModel
            The data generating model.
        single : bool
            If True, run the sigmasep algorithm once.
            If False, run the sigmasep algorithm once per possible feature (edges and confounders)
            with additional hard constraints.
        save_CIR : bool
            If True, save the conditional independence relations to file.
        alpha : float
            The significance level for the conditional independence tests for ASP.
        '''
        # Run ASP algorithm
        logger.info('Running ASP mode %d ...', mode)
        # Transform CIRs to ASP encoding
        asp = Path(out_path, 'asp.pl')
        transform_tests_file_to_asp_file(data, n_obs, asp, use_CIT, alpha=alpha)
        times = run(out_path, asp_path, asp, d=n_obs, seps=seps, mode=mode)

        return times

    @staticmethod
    def create_data(model: LinearModel, n_samples: int, exp_list: List[List[int]]) -> Dict[Tuple[int, ...], np.ndarray]:
        '''
        Samples data from model.


        Return:
        ------
        Dictionary of data. The keys are tuples indicating the performed interventions and the values are the samples.

        Parameters:
        ----------
        model: LinearModel
            Model to sample.
        n_samples: int
            Number of samples per intervention.
        exp_list: List[List[int]]
            List of list of indexes of variables to intervene on.
            Example: [[], [0], [1], [2]] stands for data from passive state,
            and three single variable interventions on nodes 0, 1 and 2 respectively.
        '''
        # complete_interventions = model.get_complete_interventions()
        # experiment_interventions = complete_interventions[:experiment + 1, :]
        exp_list_binary = index_to_binary(exp_list, model.graph.n_obs)
        data = {tuple(interventions): model.sample(n_samples=n_samples, interventions=interventions)
                for interventions in exp_list_binary}
        return data

    @staticmethod
    def save_model(model: LinearModel, path: Path) -> None:
        '''
        Saves adjacency matrix, confounder matrix and edge matrix of a linear model to file.


        Return:
        ------
        None

        Parameters:
        ----------
        model: LinearModel
            Model to save.
        path: Path
            Directory were to save data.
        '''
        path.mkdir(exist_ok=True)
        np.savetxt(str(Path(path, 'adj_true.csv').absolute()), model.graph.adjacency, delimiter=',', fmt='%d')
        np.savetxt(str(Path(path, 'confs_true.csv').absolute()), model.graph.confs, delimiter=',', fmt='%d')
        np.savetxt(str(Path(path, 'edges_true.csv').absolute()), model.graph.adjacency_observed, delimiter=',', fmt='%d')

    @staticmethod
    def save_llc(data: Dict, path: Path) -> Path:
        '''
        Saves data to pickle file.


        Return:
        ------
        The path were the data is saved.

        Parameters:
        ---------
        data: any
            Data to be pickled.
        path: Path
            Directory were to save data.
        '''
        p = Path(path, 'data.pkl')
        with p.open('wb') as f:
            pickle.dump(data, f)
        return p
