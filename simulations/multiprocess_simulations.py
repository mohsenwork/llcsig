from multiprocessing import Pool
from .app import Simulations
from itertools import product
import yaml


class Multiprocessing:

    @staticmethod
    def task(i, config):
        # One task consists of simulating one model
        Simulations.run(**config, models=[i])


    @staticmethod
    def run(*args):
        # parse configuration file
        configs = yaml.safe_load(open(f"config/{args[0]}"))

        # parse number of models
        n_models = int(args[1])

        # create a pool of workers that automatically uses all CPU cores
        with Pool() as pool:

            # issue one task for each call to the function, tasks are divided upon workers
            pool.starmap(Multiprocessing.task, product(range(100, 100 + n_models), configs))