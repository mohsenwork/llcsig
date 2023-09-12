from copy import deepcopy
from itertools import product
import numpy as np
from simulations.data_generation.directed_graph import DirectedGraph


def random_coefficients_stable(adjacency_matrix: np.ndarray, coefficient_low: float = 0.1,
                               coefficient_high: float = 1.1) -> np.ndarray:
    '''
    Creates a stable coefficient matrix for a given adjacency matrix.



    Return:
    ------
    A coefficient matrix, the same size as the adjacency matrix.

    Matrix construction:
    --------------------
    See random_coefficients(). The coefficient matrix is sampled until a stable one is found.

    Parameters:
    ----------
    adjacency : np.ndarray
        A square adjacency matrix.
    coefficient_low : float
        A positive number for the lower absolute threshold of the edge coefficients.
    coefficient_high : float
        A positive number for the upeer absolute threshold of the edge coefficients.

    '''
    while True:
        coefficient_matrix = random_coefficients(adjacency_matrix, coefficient_low, coefficient_high)

        if is_stable(coefficient_matrix):
            break
    return coefficient_matrix


def is_stable(coefficient_matrix: np.ndarray, eigen_theshold: float = 0.9) -> np.ndarray:
    '''
    Tests the asymptotic stability of the coeffiecient matrix.
    The coeffiecient matrix B is stable iff for every possible intervention indicator U (including the passive observational state)
    the eigenvalues of U*B are in a unit circle.
    If B.shape > 7, only test U=I for computational reasons.


    Return:
    ------
    True if the coefficient matrix is stable.

    Parameters:
    ----------
    eigen_threshold : float
        The threshold to test the absolute values of the eigenvalues against (default is 1.0).

    '''
    n = coefficient_matrix.shape[0]
    if n > 7:
        all_interventions = [[1]*n]
    else:
        all_interventions = list(product([0, 1], repeat=n))

    for intervention in all_interventions:

        U = np.diag(intervention)
        UB = U @ coefficient_matrix
        eigen_vals, _ = np.linalg.eig(UB)

        if np.amax(np.abs(eigen_vals)) > eigen_theshold:
            return False

    return True


def random_coefficients(adjacency_matrix: np.ndarray, coefficient_low: float = 0.1,
                        coefficient_high: float = 1.1) -> np.ndarray:
    '''
    Creates a (potentially unstable) coefficient matrix for a given adjacency matrix.


    Return:
    ------
    A coefficient matrix, the same size as the adjacency matrix.

    Matrix construction:
    --------------------
    The edge coefficients are sampled from a uniform distrubition in the
    intervals [-coefficient_high, -coefficient_low] and [+coefficient_low, +coefficient_high].

    Parameters:
    ----------
    adjacency : np.ndarray
        A square adjacency matrix.
    coefficient_low : float
        A positive number for the lower absolute threshold of the edge coefficients.
    coefficient_high : float
        A positive number for the upeer absolute threshold of the edge coefficients.

    '''
    n = adjacency_matrix.shape[0]
    coefficients_pos = np.random.uniform(low=coefficient_low, high=coefficient_high, size=(n, n))
    mask = np.random.choice([-1, 1], size=(n, n))
    coefficients = np.multiply(coefficients_pos, mask)
    coefficient_matrix = np.multiply(adjacency_matrix, coefficients)

    return coefficient_matrix


class LinearModel:

    @property
    def graph(self) -> DirectedGraph:
        # Directed graph corresponding to the model
        return self._graph

    @property
    def coefficient_matrix(self) -> np.ndarray:
        # Coefficient matrix of latent and observed
        return deepcopy(self._coefficient_matrix)

    @property
    def disturbance_std(self) -> np.ndarray:
        # Matrix of disturbance standard diviation
        return deepcopy(self._disturbance_std)

    @property
    def disturbance_mean(self) -> np.ndarray:
        # Matrix of disturbance means
        return deepcopy(self._disturbance_mean)

    def __init__(self, n_obs: int, n_conf: int, cyclic: bool = False, lower_bound: int = 0,
                 upper_bound: int = 2, indegree: bool = True, std: float = 1):
        '''
        Generates a random linear model.


        Model construction:
        --------------------
        The underlying graph is constructed using /data_generation/directed-graph.py.
        The coefficient matrix is sampled randomly until a stable model is found.
        The disturbances have zero mean.

        Parameters:-
        ----------
        n_obs : int
            The number of observed variables.
        n_conf : int
            The number of latent variable confounders. Each latent variable has exactly observed two children.
        cyclic : bool
            True if underlying graph should be cyclic.
        lower_bound : int
            The lower limit for the number of outgoing or ingoing edges, depending on indegree.
        upper_bound : int
            The upper limit for the number of outgoing or ingoing edges, depending on indegree.
        indegree : bool
            If True, then the lower and upper bound parameters specify the indegree, else outdegree.
        std : float
            Standard diviation factor.
            The standard diviations of n error variables is drawn from np.abs(np.random.normal(size=n)*0.5)+1 then multiplied by std.
            If you want the standard diviations to be small, set this parameter to a small number. 
        '''
        self._graph = DirectedGraph(n_obs, n_conf, cyclic, lower_bound, upper_bound, indegree)
        self._coefficient_matrix = random_coefficients_stable(self.graph.adjacency)
        self._disturbance_std = (np.abs(np.random.normal(size=n_obs + n_conf) * 0.5) + 1) * std
        self._disturbance_mean = np.zeros(n_obs + n_conf)

    def sample(self, n_samples: int, interventions: np.ndarray) -> np.ndarray:
        '''
        Generates data from a linear model with the given coefficient matrix.


        Return:
        ------
        A (n_obs x n_samples) matrix of samples for the observed variables.

        Parameters:
        ----------
        n_samples : int
            The number of samples.
        n_obs : int
            The number of observed variables.
        interventions : list or numpy.ndarray
            A binary vector indicating whether the observed variable should be intereved upon.
            An entry of 1 indicates intervention, else 0.

        '''
        indicator_observed = self._get_indicator_observed(interventions)
        error_matrix = self._get_error_matrix(interventions, n_samples)
        is_observed = np.array([0] * self.graph.n_conf + [1] * self.graph.n_obs, dtype=bool)
        X = np.linalg.inv(np.eye(self.graph.n_obs + self.graph.n_conf) - (indicator_observed @ self.coefficient_matrix)) @ error_matrix

        return X[is_observed]

    def get_random_interventions(self, n_intervened: int) -> np.ndarray:
        '''
        Creates an indicator array of random interventions on n_intervened variables out of n_obs variables.


        Return:
        ------
        A binary vector of length n_obs indicating whether the observed variable should be intereved upon.
        An entry of 1 indicates intervention, else 0.

        Parameters:
        ----------
        n_intervened : int
            The number of interventions.

        '''
        interventions = np.zeros(self.graph.n_obs)
        interventions[:n_intervened] = 1
        np.random.shuffle(interventions)

        return interventions.astype(int)

    def get_complete_interventions(self) -> np.ndarray:
        '''
        Creates indicator arrays of all single variable interventions and zero variable intervention.


        Return:
        ------
        An array consisting of n_obs indicator arrays.
        An indicator array is a binary vector of length n_obs indicating whether the observed variable
        should be intereved upon. An entry of 1 indicates intervention, else 0. The first row is for the
        passive observational state (all 0s), the second for row for first observed variable and so on...

        '''
        return np.vstack([np.zeros(self.graph.n_obs), np.eye(self.graph.n_obs)]).astype(int)

    def _get_error_matrix(self, interventions: np.ndarray, n_samples: int) -> np.ndarray:
        '''
        Creates samples of error matrix for data generation.
        Each row contains the errors of one variable (confounder then observed), which are independent.
        Intervened variables have zero mean and unit variance.


        Return:
        ------
        A (n_conf + n_obs x n_samples) matrix.

        Parameters:
        ----------
        interventions : np.ndarray
            A binary vector indicating whether the observed variable should be intereved upon.
            See LinearModel.get_complete_interventions() or LinearModel.get_random_interventions().
        n_samples : int
            The number of samples.
        '''
        is_intervened = np.array(np.concatenate([[0] * self.graph.n_conf, interventions]), dtype=bool)
        cov = np.diag(np.power(self.disturbance_std, 2))
        cov[is_intervened, is_intervened] = 1
        mean = self.disturbance_mean
        mean[is_intervened] = 0

        return np.random.multivariate_normal(mean=mean, cov=cov, size=n_samples).T

    def _get_indicator_observed(self, interventions: np.ndarray) -> np.ndarray:
        '''
        Creates diagonal indicator matrix U where U[i,i] = 1 iff variable i is passively observed.


        Return:
        ------
        A (n_obs+n_conf x n_obs+n_conf) square matrix.

        Parameters:
        ----------
        interventions : np.ndarray
            A binary vector indicating whether the observed variable should be intereved upon.
            See LinearModel.get_complete_interventions() or LinearModel.get_random_interventions().
        '''
        indicator_intervened = np.diag(np.concatenate([[0] * self.graph.n_conf, interventions]))
        return np.eye(self.graph.n_obs + self.graph.n_conf) - indicator_intervened

    def _get_disturbance_covariance_matrix(self) -> np.ndarray:
        '''
        Calculates the covariance of the errors/disturbances Σe (as described in LLC paper), where Σe(x_i, x_j) = cov(e_i, e_j), e_i being the
        disturbance term corresponding to x_i and e_j corresponding to x_j. A non-zero covariance means confounding. (Lemma 7)

        Return:
        ------
        A (n_obs x n_obs) matrix.
        '''
        B_obs = self.coefficient_matrix[self.graph.n_conf: ,self.graph.n_conf:]
        I_obs = np.eye(self.graph.n_obs)
        T = np.linalg.inv(np.eye(self.graph.n_obs + self.graph.n_conf) - self.coefficient_matrix)
        E = np.diag(np.power(self.disturbance_std, 2))
        Cx = (T @ E @ T.T)[self.graph.n_conf: ,self.graph.n_conf:]
        E_obs = (I_obs - B_obs) @ Cx @ (I_obs - B_obs).T
        return E_obs


    def _get_covariance_observed_infinite(self, interventions: np.ndarray) -> np.ndarray:
        '''
        Calculates covariance of observed variables from B and Σe (as described in LLC paper).
        For an intervention E = (J, U) the covariance is (I-UB)^-1(J+UΣeU)(I-UB)^-T.
        Where A^-1 is the inverse of A and A^-T is the inverse of the transpose.
        Marginalise latent variables by only returning the covariance for observed nodes.

        Return:
        ------
        A (n_obs x n_obs) matrix.

        Parameters:
        ----------
        interventions : np.ndarray
            A binary vector indicating whether the observed variable should be intereved upon.
            1 means intervention.
            See LinearModel.get_complete_interventions() or LinearModel.get_random_interventions().
        '''
        U = self._get_indicator_observed(interventions)
        I = np.eye(self.graph.n_obs + self.graph.n_conf)
        J = I - U
        B = self.coefficient_matrix
        T = np.linalg.inv(I - U @ B)
        E = np.diag(np.power(self.disturbance_std, 2))
        C = T @ (J + U @ E @ U) @ T.T
        return C[self.graph.n_conf: ,self.graph.n_conf:]
