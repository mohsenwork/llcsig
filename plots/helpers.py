import re
import sys
import os
import pickle
import numpy as np
from pathlib import Path
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

def roc_pr_stats(reader, filepath, edge_type='edge', sample_size=10000, experiment='5', mode=2, faithful='nf', sep='s'):
    '''
    Generates stats about the ROC curve.
    '''
    scores = np.empty(shape=0)
    labels = np.empty(shape=0)
    models = sorted([ item for item in os.listdir(filepath) if re.match('model_*', item) ], key=lambda s: int(s[6:]))

    for model in models:

        path = Path(filepath, f'{model}', f'samplesize_{sample_size}', f'experiment_{experiment}')
        A_pred = reader(path, edge_type, mode, faithful, sep)
        path_true = Path(filepath, f'{model}', f'{edge_type}s_true.csv')
        A_true = np.loadtxt(path_true, delimiter=',')
        d = A_true.shape[0]

        for row in range(d):
            for col in range(d):
                if row != col:
                    labels = np.append(labels,int(A_true[row,col]))
                    scores = np.append(scores,A_pred[row,col])

    fpr, tpr, _ = roc_curve(y_true=labels, y_score=scores, pos_label=[1])
    roc_auc = auc(fpr, tpr)
    prec1, rec1, _ = precision_recall_curve(y_true=labels, probas_pred=scores, pos_label=[1])
    prec0, rec0, _ = precision_recall_curve(y_true=labels, probas_pred=-scores, pos_label=[0])
    avg_prec1 = average_precision_score(labels, scores)
    avg_prec0 = average_precision_score(labels, -scores)
    return (fpr,tpr,roc_auc,prec0,rec0,prec1,rec1,avg_prec0,avg_prec1)


def read_llc(path, edge_type='edge', mode=2, faithful='nf', sep='s'):
    '''
    Helper for loading llc files for roc_pr_stats.
    '''
    path_pred = Path(path, f'llc_{faithful}_{mode}', 'llc_out.pkl')

    with open(path_pred, 'rb') as f:
        data = pickle.load(f)

    if mode == 2:
        if edge_type == "edge":
            return np.abs(data['Bz'])
        else:
            return np.abs(np.nan_to_num(data['Cez'], nan=0))
    else: 
        if edge_type == "edge":
            return np.abs(data['B'])
        else:
            return np.abs(data['Ce'])
        

def read_asp(path, edge_type, mode=2, faithful='nf', sep='s'):
    '''
    Helper for loading asp files for roc_pr_stats.
    '''
    path_pred = Path(path, f'asp_{mode}', f'{edge_type}s_score_{sep}_sep.csv')
    return np.loadtxt(path_pred, delimiter=',')


def mask(A, threshold=5):
    '''
    Returns a copy of A where elements below/equal to the threshold are set to 0 and above to 1.
    '''
    _A = np.copy(A)
    _A[A <= threshold] = 0
    _A[A > threshold] = 1
    return _A


def identification_stats(A_true_, A_pred_, type='edge'):
    '''
    Takes the true and predicted feature matrix as input and computes a confusion matrix of the form
                                          Actual
                        positive               negative    
              positive  true positive (TP)     false positive (FP)
    Predicted unknown   positive unknwon (PU)  negative unknown (NU)
              negative  false negative (FN)    true negative (TN)
    
    where positive = feature is present, negative = feature is not present
    returns an array with [TP, FP, PU, NU, FN, TN]
    '''    
    # Set diagonals to 0
    A_true = np.copy(A_true_)
    np.fill_diagonal(A_true, 0)
    A_pred = np.copy(A_pred_)
    np.fill_diagonal(A_pred, 0)
 
    ## true 
    TP = np.count_nonzero(np.logical_and(A_pred == 1, A_true == 1))
    tmp = np.logical_and(A_pred == 0, A_true == 0)
    np.fill_diagonal(tmp, False)
    TN = np.count_nonzero(tmp)

    ## unknown
    PU = np.count_nonzero(np.logical_and(A_pred == 0.5, A_true == 1))
    NU = np.count_nonzero(np.logical_and(A_pred == 0.5, A_true == 0))

    ## false
    FN = np.count_nonzero(np.logical_and(A_pred == 0, A_true == 1))
    FP = np.count_nonzero(np.logical_and(A_pred == 1, A_true == 0))

    # divide stats about confounders by 2
    if type == 'conf':
        TP /= 2; TN /= 2; PU /= 2; NU /= 2; FN /= 2; FP /= 2;  

    return [TP, FP, PU, NU, FN, TN]

def accuracy(TP, FP, FN, TN):
    '''
    Calculates the accuracy from true positive, false positive, false negative and true negative values.
    '''
    return (TP+TN) / (TP+TN+FP+FN)

def f1_score(TP, FP, FN, TN):
    '''
    Calculates the F1 score from true positive, false positive, false negative and true negative values.
    '''
    return (2*TP) / (2*TP + FP + FN)

def get_true(filepath):
    '''
    Returns true matrices of created models.
    1 means feature is present
    0 means feature is absent
    Element i of the output describes model i. The first element of that pair is the 
    matrix for the edges, the second for the confounders.
    '''
    models = sorted([ item for item in os.listdir(filepath) if re.match('model_*', item) ], key=lambda s: int(s[6:]))
    out = []
 
    for model in models:
        path_edge = Path(filepath, f'{model}', 'edges_true.csv')
        path_conf = Path(filepath, f'{model}', 'confs_true.csv')
        true_edge = np.loadtxt(path_edge, delimiter=',')
        true_conf = np.loadtxt(path_conf, delimiter=',')

        out.append((true_edge, true_conf))

    return out