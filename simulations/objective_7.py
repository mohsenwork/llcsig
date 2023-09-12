import os, re
from pathlib import Path
from .app import Simulations
from .objective import Objective
from plots.helpers import accuracy, identification_stats, read_asp, mask
import numpy as np
from sklearn.metrics import auc, roc_curve
from itertools import product


class Objective7(Objective):

    def __call__(self, trial):
        """
        Evaluates a set of parameters and returns a score.
        The hyperparameter is lambda ranging between 0 and 1.
        
        Parameters:
            trial (Trial): An instance of the Trial class from Optuna.
            
        Returns:
            float: The score for the set of parameters.
        """
        alpha = trial.suggest_float("alpha", 0.005, 0.2)
        z   = trial.suggest_int("z", - 10, 50)

        for i in range(self.n_models):
            self.configs[i]['name'] = str(self.path.relative_to(self.path.parent.parent.parent))
            Simulations.run(**self.configs[i], models=[i], alpha_asp=alpha)

        # Calculate an objective value
        roc_auc = Objective7.metric(self.path, self.configs[0]['n_samples'], self.configs[0]['seps'][0], z)
        
        return roc_auc

    @staticmethod
    def metric(path, sample_sizes, sep, z):
        """
        This function calculates the total accuracy of edges and confounders from confidence 
        values over all experiments and sample sizes. 

        Parameters:
        path (path): specifies where the results of 'LLC' are saved.
        sample_sizes (List[Union[int, str]]): specifies the used sample sizes.
        sep (str: {'s', 'd'}): specifies the ASP-version.
        z (int): specifies the threshold for the confidence values.


        Returns:
        float: total roc auc of edges and confounders over all experiments and sample sizes.
        """
        TP = FP = FN = TN = 0;    
        edge_types = ['edge', 'conf']
        experiments = ['1', '2']
        models = sorted([ item for item in os.listdir(path) if re.match('model_*', item) ], key=lambda s: int(s[6:]))

        for model in models:
            for edge_type, experiment, sample_size in product(edge_types, experiments, sample_sizes):

                # read confidence values
                dir_path = Path(path, f'{model}', f'samplesize_{sample_size}', f'experiment_{experiment}')
                A_pred = read_asp(dir_path, edge_type, mode=3, sep=sep)
                path_true = Path(path, f'{model}', f'{edge_type}s_true.csv')
                A_true = np.loadtxt(path_true, delimiter=',')

                # apply threshold
                A_pred_masked = mask(A_pred, threshold=z)

                # identification stats
                tp, fp, pu, nu, fn, tn = identification_stats(A_true, A_pred_masked, type='edge')
                assert pu == 0
                assert nu == 0

                # aggregate stats
                TP += tp
                FP += fp
                FN += fn
                TN += tn
        
        acc = accuracy(TP, FP, FN, TN)
        return acc
        
