import os, re
from pathlib import Path
from .app import Simulations
from .objective import Objective
from plots.helpers import read_llc, identification_stats, accuracy, mask
import numpy as np
from itertools import product


class Objective2(Objective):

    def __call__(self, trial):
        """
        Evaluates a set of parameters and returns a score.
        The hyperparameter is lambda ranging between 0 and 1.
        
        Parameters:
            trial (Trial): An instance of the Trial class from Optuna.
            
        Returns:
            float: The score for the set of parameters.
        """
        reg = trial.suggest_float("lambda", 0.0, 1.0)
        z   = trial.suggest_float("z", 0.0, 20)

        for i in range(self.n_models):
            self.configs[i]['name'] = str(self.path.relative_to(self.path.parent.parent.parent))
            Simulations.run(**self.configs[i], models=[i], reg=reg)

        # Calculate an objective value
        roc_auc = Objective2.metric(self.path, self.configs[0]['n_samples'], z)
        
        return roc_auc

    @staticmethod
    def metric(path, sample_sizes, z):
        """
        This function calculates the total accuracy of edges and confounders from z-scores 
        over all experiments and sample sizes. 

        Parameters:
        path (path): specifies where the results of 'LLC' are saved.
        sample_sizes (List[Union[int, str]]): specifies the used sample sizes.
        z (float): specifies the used z-score threshold.


        Returns:
        float: total roc auc of edges and confounders over all experiments and sample sizes.
        """
        TP = FP = FN = TN = 0;
        edge_types = ['edge', 'conf']
        experiments = ['1', '2']
        models = sorted([ item for item in os.listdir(path) if re.match('model_*', item) ], key=lambda s: int(s[6:]))

        for model in models:
            for edge_type, experiment, sample_size in product(edge_types, experiments, sample_sizes):

                # calculate z-score
                dir_path = Path(path, f'{model}', f'samplesize_{sample_size}', f'experiment_{experiment}')
                A_pred = read_llc(dir_path, edge_type, mode=2, faithful='nf', sep='')
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
