import os, re
from pathlib import Path
from .app import Simulations
from .objective import Objective
from plots.helpers import roc_pr_stats, read_asp
import numpy as np
from sklearn.metrics import auc, roc_curve
from itertools import product


class Objective5(Objective):

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

        for i in range(self.n_models):
            self.configs[i]['name'] = str(self.path.relative_to(self.path.parent.parent.parent))
            Simulations.run(**self.configs[i], models=[i], alpha_asp=alpha)

        # Calculate an objective value
        roc_auc = Objective5.metric(self.path, self.configs[0]['n_samples'], self.configs[0]['seps'][0])
        
        return roc_auc

    @staticmethod
    def metric(path, sample_sizes, sep):
        """
        This function calculates the total roc auc of edges and confounders from z-scores 
        over all experiments and sample sizes. 

        Parameters:
        path (path): specifies where the results of 'LLC' are saved.
        sample_sizes (List[Union[int, str]]): specifies the used sample sizes.
        sep (str: {'s', 'd'}): specifies the ASP-version.

        Returns:
        float: total roc auc of edges and confounders over all experiments and sample sizes.
        """
        scores = np.empty(shape=0)
        labels = np.empty(shape=0)
        edge_types = ['edge', 'conf']
        experiments = ['1', '2']
        models = sorted([ item for item in os.listdir(path) if re.match('model_*', item) ], key=lambda s: int(s[6:]))

        for model in models:
            for edge_type, experiment, sample_size in product(edge_types, experiments, sample_sizes):

                dir_path = Path(path, f'{model}', f'samplesize_{sample_size}', f'experiment_{experiment}')

                A_pred = read_asp(dir_path, edge_type, mode=3, sep=sep)
                
                path_true = Path(path, f'{model}', f'{edge_type}s_true.csv')
                A_true = np.loadtxt(path_true, delimiter=',')
                d = A_true.shape[0]

                for row in range(d):
                    for col in range(d):
                        
                        if row != col:
                            labels = np.append(labels,int(A_true[row,col]))
                            scores = np.append(scores,A_pred[row,col])
        
        fpr, tpr, _ = roc_curve(y_true=labels, y_score=scores, pos_label=[1])
        roc_auc = auc(fpr, tpr)
        return roc_auc
        
