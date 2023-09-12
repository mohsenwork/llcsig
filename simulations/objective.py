from abc import ABC, abstractmethod, abstractstaticmethod


class Objective(ABC):
    """
    Specifies the objective function for an Optuna study.
    
    Attributes:
        configs (array): An array of dictionaries that specifies the parameters for the function Simulations.run().
                         There is one dictionary for each model.
        trial_num (int): The number of the trial.
        n_models (int): The number of models to be simulated per trial.
        
    Methods:
        __call__: Implements the objective function that evaluates a set of parameters and returns a score.
        metric: Calculates the metric optimized over in an optuna hyperparameter optimization study.
    """
    def __init__(self, configs, n_models, path):
        # Class initialization code
        self.configs = configs
        self.n_models = n_models
        self.path = path

    @abstractmethod
    def __call__(self, trial):
        """
        Evaluates a set of parameters and returns a score.
        
        Parameters:
            trial (Trial): An instance of the Trial class from Optuna.
            
        Returns:
            float: The score for the set of parameters.
        """
        pass

    @abstractstaticmethod
    def metric(self, *args):
        """
        Calculates the metric optimized over in an optuna hyperparameter optimization study.
        
        Parameters:
            args: arbitrary number of arguments.
        
        Returns:
            float: the metric.
        """
        pass