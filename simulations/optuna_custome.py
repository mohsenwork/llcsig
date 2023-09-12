from pathlib import Path
import os
import optuna
import yaml
from multiprocessing import Pool, cpu_count
from .objective import Objective

class OptunaCustome:

    @staticmethod
    def task(class_objective: Objective, study_name, storage_name, configs, n_models, path):
        """
        Defines and executes a single task as part of an Optuna study.

        Parameters:
            study_name (str): The name of the Optuna study.
            storage_name (str): The storage name for the Optuna study in the form of a relational database.
            configs (array): An array of dictionaries that specifies the parameters for the function Simulations.run().
            n_models (int): The number of models to be simulated per trial.
            path (path): The path where results are saved.
            class_objective (class): A class implementing the abstract Objective class.

        Returns:
            None
        """
        # load load_study
        study = optuna.load_study(study_name=study_name, storage=storage_name)
        study.optimize(class_objective(configs, n_models, path), n_trials=1)

    @staticmethod
    def run(class_objective: Objective, *args):
        # parse input arguments
        study_name = args[0]
        configs = yaml.safe_load(open(f"config/{args[1]}"))
        n_models = len(configs)
        n_trials = int(args[2])

        # study and storage names
        path = Path(Path.cwd(), 'simulations', 'data', 'optuna', study_name)
        path.mkdir(parents=True, exist_ok=True)
        storage_name = f"sqlite:///{path}/{study_name}.db"

        # create study
        optuna.create_study(study_name=study_name, storage=storage_name, sampler=optuna.samplers.TPESampler(seed=1), direction="maximize", load_if_exists=True)

        # get last completed trial
        trial_folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f)) and f.startswith("trial_")]
        trial_numbers = [int(f.split("_")[1]) for f in trial_folders]
        max_trial_number = max(trial_numbers, default=-1) + 1

        # create a pool of workers that automatically uses all CPU cores
        with Pool(cpu_count()) as pool:

            # issue one task for each call to the function, tasks are divided upon workers
            pool.starmap(OptunaCustome.task, [(class_objective, study_name, storage_name, configs, n_models, Path(path, f'trial_{trial_num}')) 
                                          for trial_num in range(max_trial_number, max_trial_number + n_trials)])