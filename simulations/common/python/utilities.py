import logging
from pathlib import Path
from typing import List, Tuple, Dict
import re
import pandas as pd
import numpy as np
import scipy.stats as stats
from simulations.resources.sigmasep.main import all_xyZ_partitions_list
from simulations.data_generation.linear_model import LinearModel
logger = logging.getLogger(__name__)


def load_globals(path: Path):
    pass


def d_sep_CIRS(model: LinearModel, I: List[List[int]]) -> pd.DataFrame:
    '''
    Computes conditional independence relations according to the d-seperation criterion.
    For a particular set of intervened variables J, each partition x,y,Z of the variables V
    is tested for d-seperation.
    If nodes x and y are d-seperated by the conditioning set Z in the graph with interventions on J,
    then the associated row in the output dataframe denotes True in the 'sep' column.

    Return:
    ------
    A Dataframe with d-seperation statements.

    Parameters:
    ----------
    model : LinearModel
        A linear model.
    I : List[List[int]]
        A list of sublists of indexes.
        Each sublist contains the indexes of intervened variables (starting at 0).
    '''
    nodes = list(range(model.graph.n_obs))
    U = np.zeros((model.graph.n_obs, model.graph.n_obs))
    B = model.graph.confs
    D = model.graph.adjacency_observed

    data = [d_sep_row(x[0], y[0], Z, do_targets, D, B, U)
            for do_targets in I
            for x, y, Z in all_xyZ_partitions_list(nodes)]

    return pd.DataFrame(data)


def directed_reachable(x, y, C, J, D, B, U=None):
    # D-separation for cyclic mixed graphs code from paper:
    # "Constraint-Based Causal Discovery: Conflict Resolution with Answer Set Programming"
    # By Antti Hyttinen, Frederick Eberhardt, Matti Järvisalo
    # The code below is a translation of the R code provided by the authors.

    # Implements a d-separation oracle.
    # x,y    variables considered
    # C      vector of variable indexes in the conditioning set
    # J      vector of variable indexes in the intervention set  
    # D      binary matrix of directed edges, G[i,j] = 1 means j->i is present
    # B      symmetric binary matrix of bidirected edges <-> , for latent confounding
    # U      symmetric binary matrix of undirected edges --- , should be empty at beginning

    # Note that this is not currently the most efficient implementation.
    # - One could use matrix operations
    # - One could calculate more of the relations all at once.

    n = D.shape[0]

    D[J] = 0
    B[J] = 0
    B[:, J] = 0

    ###########################################################################

    HH = B  # symmetric head head paths

    if U is None:  # there should generally be no head-head paths in the beginning 
        TT = np.zeros((n, n))  # this is for tail tail paths
    else:
        TT = U
    TH = D  # notice here TH[x,y] = 1 iff path y->x

    # TH and HH self paths do not affect the d-connectedness so they can be ignored
    np.fill_diagonal(TH, -1)
    np.fill_diagonal(HH, -1)

    for node in range(n):    # doing either conditioning or marginalizing to all variables

        if node == x or node == y:
            continue  # skip variables that are in the final graph

        # gather up all different kinds of parents, or connected nodes
        thpa = np.where(TH[node, ] == 1)[0]
        htpa = np.where(TH[:, node] == 1)[0]  # aka children
        hhpa = np.where(HH[node, ] == 1)[0]
        ttpa = np.where(TT[node, ] == 1)[0]

        if node not in C:
            # the marginalization operation is more difficult
            # i-->node-->j
            for i in thpa:
                for j in htpa:
                    TH[j, i] = 1

            # i--> node --- j
            for i in thpa:
                for j in ttpa:
                    TT[i, j] = TT[j, i] = 1

            # i --> node <-- is not ok
            # i --> node <-> is not ok
            #####################################

            # i <-> node --> j
            for i in hhpa:
                for j in htpa:
                    HH[j, i] = HH[i, j] = 1

            # i <-> node --- j
            for i in hhpa:
                for j in ttpa:
                    TH[i, j] = 1

            # i <-> node <-- is not ok
            # i <-> node <-> is not ok
            ######################################

            # i --- node --- j
            # tail tail parents connected
            for i in ttpa:  # connects node to itself as well so tt-self cycle is inherited
                for j in ttpa:
                    TT[i, j] = TT[j, i] = 1
            # i --- node --> j
            for i in ttpa:  # connects node to itself as well so tt-self cycle is inherited
                for j in htpa:
                    TH[j, i] = 1
            # i --- node <-> j done already
            # i --- node <-- j done already
            ##############################################

            # i <-- node --> j
            for i in htpa:  # connects node to itself as well so tt-self cycle is inherited
                for j in htpa:
                    HH[i, j] = HH[j, i] = 1
            # i <-- node <-> j done already
            # for i in htpa:  #connects node to itself as well so tt-self cycle is inherited
            #  for j in hhpa:
            #    HH[i,j]<-HH[j,i]<-1  
            # }    
            # i --- node <-> j done already
            # i --- node <-- j done already

        if node in C or TT[node, node] == 1:
            # notice the simplicity here!
            # an unconditioned node with a selfloop actually allows through all traffic!!!
            # only three options
            # i--> node <--j
            for i in thpa:
                for j in thpa:  # notice that this connects that parents to them selves as well
                    TT[i, j] = TT[j, i] = 1
            # i<-> node <->j
            # hh parents need to be connected by head head nodes
            for i in hhpa:
                for j in hhpa:
                    HH[i, j] = HH[j, i] = 1
            # i<-> node <--j
            # connecting hh parent to th parent
            for i in hhpa:
                for j in thpa:
                    TH[i, j] = 1

        # only tailtail cycles are relevant to d-connection
        np.fill_diagonal(TH, -1)
        np.fill_diagonal(HH, -1)

        # now take the node away
        TH[node] = TH[:, node] = HH[:, node] = HH[node, :] = TT[:, node] = TT[node, :] = -1

    # the nodes are connected if any of the paths is present in the end, where
    # all variables 
    return (HH[x, y] == 1 or TH[x, y] == 1 or TH[y, x] == 1 or HH[y, x] == 1 or TT[y, x] == 1)


def d_sep_row(x: int, y: int, C: List[int], J: List[int], D: np.array, B: np.array, U: np.array) -> Dict:
    '''
    Helper; computes for CIR in dictionary form.
    '''
    sep = d_sep(x, y, C, J, D, B, U)
    return {'X': f'[{x}]',
            'Y': f'[{y}]',
            'do_targets': f'[{";".join([str(s) for s in J])}]',
            'Z': f'[{";".join([str(s) for s in C])}]',
            'sep': sep}


def d_sep(x: int, y: int, C: List[int], J: List[int], D: np.array, B: np.array, U: np.array) -> bool:
    '''
    Test d-separation of x and y conditioned on C and with interventions J in the model specified by D, B and U.

    Return:
    ------
    True, if the nodes are seperated.

    Parameters:
    ----------
    x : int
        Index of x.
    y : int
        Index of y.
    C : List[int]
        Index list of conditioning set.
    J : List[int]
        Index list of intervention set.
    D : np.array
        A binary matrix of directed edges, D[i,j] = 1 iff j->i is present.
    B : np.array
        A symmetric binary matrix of bidirected edges, D[i,j] = 1 iff j and i are confounded.
    U : np.array
        A symmetric binary matrix of undirected edges, This is usually be empty at beginning.
    '''
    _C = np.copy(C).astype(int)
    _J = np.copy(J).astype(int)
    _D = np.copy(D).astype(int)
    _B = np.copy(B).astype(int)
    _U = np.copy(U).astype(int)
    out = directed_reachable(x, y, _C, _J, _D, _B, _U)
    return not out


def index_to_binary(exp_list: List[List[int]], n: int) -> List[List[int]]:
    '''
    Translates intervention notations from indexes to binary.
    Example  [[], [0], [1]] to [[0,0,0,0,0], [1,0,0,0,0],[0,1,0,0,0]].
    '''
    out = []
    for exp in exp_list:
        l = [0]*n
        for e in exp:
            l[e] = 1
        out.append(l)
    return out


def pcor_dataframe(dic):
    '''
    Computes partial correlation coefficients from the covariance matrix.
    For a particular set of intervened variables J (keys of dic), for each partition x,y,Z of the variables V
    the partial correlation coefficients are computed from the given covariance matrix (values of dic).
    Since computation of partial correlation coefficients from a covariance matrix requries it to be positive definite - 
    which is not always the case - a value np.nan is returned for such cases.

    Return:
    ------
    A Dataframe with partial correlation coefficients.

    Parameters:
    ----------
    dic:
        A key is a binary representation of an intervention. (E.g. [0,0,1,0])
        A value is the corresponding covariance matrix.
    '''
    rows = []

    for key in dic.keys():
        cov = dic[key]
        pos_def = is_symm_pd(cov)
        nodes = list(range(len(key)))

        for x, y, C in all_xyZ_partitions_list(nodes):
            J = np.nonzero(key)[0]
            if not pos_def:
                r = np.nan
            else:
                r = pcor(x[0], y[0], C, cov)

            rows += [{'X': f'[{x[0]}]',
                    'Y': f'[{y[0]}]',
                    'do_targets': f'[{";".join([str(s) for s in J])}]',
                    'Z': f'[{";".join([str(s) for s in C])}]',
                    'r': r}]

    return pd.DataFrame(rows)


def pcor(x: int, y: int, C: List[int], cov: np.ndarray) -> float:
    '''
    Computes the partial correlation between x and y given C from the covariance matrix.
    Based on inversion of covariance matrix.
    Python implementation of: https://github.com/cran/ggm/blob/master/R/functions.R (Lines 2784-2848)

    Return:
    ------
    A scalar.

    Parameters:
    ----------
    x: int,
        First index of variable for which the correlation is computed.
    y: int,
        Second index of variable for which the correlation is computed.
    C: List[int],
        List of indexes for the conditioning set.
    cov: np.ndarray
        Covariance matrix that is positive definite.
    '''
    u = [x, y] + C
    P = np.linalg.inv(cov[np.ix_(u, u)])
    p_cor = -P[0,1] / np.sqrt(P[0,0] * P[1,1])
    return p_cor


def pcor_test(r: int, q: int, n: int) -> Tuple[float, int, float]:
    '''
    Computes t statistic and (two-sided) p value from partial correlation coefficient.


    Return:
    ------
    A scalar

    Parameters:
    ----------
    r: int, 
        The partial correlation coefficient computed by pcor().
    q: int,
        Number of variables in the conditioning set.
    n: int
        The sample size.
    '''
    dof  = n - 2 - q # degrees of freedom
    t_stat = r * np.sqrt(dof) / np.sqrt(1 - r*r)
    p_val = 2 * stats.t.cdf(-abs(t_stat), dof) 
    return t_stat, dof, p_val


def independent(r: int, cond: int = 0, n: int = 1e15, alpha: float = 0.05) -> bool:
    '''
    Returns True if the p-value of the partial correlation of x and y given C is above the significance value alpha.
    If the partial correlation was computed form an oracle, you can use the default value for n = 1e15 for the p-value
    calculation.

    Return:
    ------
    True if independent.

    Parameters:
    ----------
    r: float,
        The correlation coefficient.
    cond: int,
        Length of the conditioning set. Used to calculated degrees of freedom.
        If the covariance matrix stemms from an oracle, use default value.
    n: int,
        Number of samples with covariance matrix was generated from.
        If the covariance matrix stemms from an oracle, use default value.
    alpha: float,
        Threshold for dependence.
    '''
    p_val = pcor_test(r, cond, n)[2]
    return p_val > alpha


def is_symm_pd(A: np.ndarray) -> bool:
    '''
    Checks if matrix is symmetric positive definite using choleskly decomposition.

    Returns:
    -------
    True, if matrix is symmetric and positive definite.

    Parameters:
    ----------
    A: np.ndarray
        A matrix.
    '''
    if np.allclose(A, A.T):
        try:
            np.linalg.cholesky(A)
            return True
        except np.linalg.LinAlgError:
            return False
    else:
        return False
    

def d_sep_from_clingo(path: Path) -> pd.DataFrame:
    '''
    Calculates d-seperation statements (with empty intervention set) from a clingo output graph.
    If nodes x and y are d-seperated by the conditioning set Z,
    then the associated row in the output dataframe denotes True in the 'sep' column.

    Return:
    ------
    A Dataframe with d-seperation statements.

    Parameters:
    ----------
    path : Path
        Path to the ASP solver output file in .txt format.
    '''
    # parse clingo output
    D, B = parse_clingo(path)
    n = D.shape[0]  # pylint: disable=E1136  # pylint/issues/3139
    U = pd.DataFrame(np.zeros((n, n)))

    # calc CIRs
    nodes = list(range(n))
    data = [d_sep_row(x[0], y[0], Z, [], D, B, U)
            for x, y, Z in all_xyZ_partitions_list(nodes)]

    return pd.DataFrame(data)


def parse_clingo(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Returns highest scoring graph (lower optimization cost) from clingo output.
    The ASP solver (clingo) outputs a sequence of answer sets.
    Each answer set contains predicates representing a solution graph along with an optimization cost.
    The optimization cost is the value sufferd by the cost function; lower is better.
    The following predicates are filtered to construct the solution graph:
        edge(X,Y), which represents an edge x → y
        conf(X,Y), which represents a confounder between x and y
        nodes(i), which represents a node i in a graph (i >= 0)

    Return:
    ------
    A tuple consisting of a adjacency matrix D and confounding matrix B.
    D[i, j] = 1 iff there is an edge from j into i.
    B[i, j] = B[j, i] = 1 iff i and j are confounded.

    Parameters:
    ----------
    path : Path
        Path to the ASP solver output file in .txt format.
    '''
    # read ASP output
    with open(path) as f:
        lines = f.readlines()

    # find answer sets along with optimization cost
    opt_indexes = [(i + 1, int(re.search('[0-9]+', lines[i + 2]).group()))
                   for i, line in enumerate(lines) if 'Answer:' in line]

    # find answer set with lowest optimization cost
    answer_index, _ = sorted(opt_indexes, key=lambda x: x[1])[0]
    answer_set = lines[answer_index]

    # filter for edges and confs predicates
    edge_predicates = re.findall(r'edge\([0-9]+,[0-9]+\)', answer_set)
    conf_predicates = re.findall(r'conf\([0-9]+,[0-9]+\)', answer_set)

    # filter for nodes predicates
    nodes = re.findall(r'node\([0-9]+\)', answer_set)
    n_obs = int(re.findall(r'[0-9]+', nodes[-1])[0])

    # create edge matrix
    D = np.zeros((n_obs + 1, n_obs + 1))
    edges = [(int(a), int(b)) for p in edge_predicates for a, b in [re.findall('[0-9]+', p)]]
    D[[b for _, b in edges], [a for a, _ in edges]] = 1

    # create conf matrix
    B = np.zeros((n_obs + 1, n_obs + 1))
    confs = [(int(a), int(b)) for p in conf_predicates for a, b in [re.findall('[0-9]+', p)]]
    B[[b for _, b in confs], [a for a, _ in confs]] = 1
    B[[a for a, _ in confs], [b for _, b in confs]] = 1

    return (D, B)
