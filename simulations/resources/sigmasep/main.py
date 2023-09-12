# pylint: skip-file
import os
from os.path import basename, dirname, splitext
from pathlib import Path
import re
from typing import Dict, List, Tuple
import numpy as np
import scipy.stats
import pandas as pd
from numpy.random import exponential
from numpy.random import choice
from numpy.random import binomial
from numpy.random import normal
from numpy.random import uniform
from numpy.linalg import norm
from itertools import product
import subprocess
from datetime import datetime
from timeit import default_timer as timer


# This file has the same logic as the original code, but allows for running clingo once and for every feature.
# Note that no changes where made to the underlying algorithm.
# Fixed outdated code:
#   Converted 4 spaces to tabs because function run_all_CIT_and_save() in 'mSCM.py' had an indentation error.
#   In function normalranktransform(), changed line
#       S = df.rank(method='first').as_matrix()/(n+1) to 
#       S = df.rank(method='first').to_numpy()/(n+1)
#   In function pcor(), added lines
#		ex = ex.reshape((ex.shape[0],))
#		ey = ey.reshape((ey.shape[0],))
#		due to compilation error in line:
#		scipy.stats.pearsonr(ex, ey)
#   because of error 'DataFrame' object has no attribute 'as_matrix' (outdated function).

def run(outdir: Path, aspdir: Path, asp_input_file: Path, d: int, seps: List[str], mode: int) -> List[float]:
	'''
	Runs the ASP algorithm. There are three different modes.


	Return :
	----------------------------------
	The output of clingo is written to a file. Returns a dictionary of runtimes; keys seps and values the runtime.

	Parameters:
	----------
	outdir : Path
		Path where the output is saved.
	aspdir : Path
		Path where the asp encoding lives.
	asp_input_file: Path
		Path for the conditional independence relations encoding file.
	d: int
		Number of observed variables
	seps : array
		Which criterion(s) to use.
		's' corresponds to sigma-separantion encoding in sigma_hej_cyclic.pl, 
		'd' corresponds to d-separantion encoding for cyclic case in hej_cyclic.pl,
		'a' corresponds to d-separantion encoding for acyclic case in hej_acyclic.pl,
	mode : int
		Which mode to run clingo in.
		mode 1: Run clingo to find one optimal answer set.
		mode 2: Run clingo once to compute the union all optimal answer sets and once to compute their intersection.
		mode 3: Run clingo once per possible features with additional hard constraints.
	''' 
	times = []

	for sep in seps:
		start = timer()

		if mode == 1:
			# create command
			out_file = Path(outdir, f'results_{sep}.txt')
			flags = '-W no-atom-undefined --quiet=1'
			show_file = Path(aspdir, 'show.pl')
			command = create_command(aspdir, [show_file, asp_input_file], out_file, sep,flags)
			# run command
			subprocess.call(command, shell=True)

		elif mode == 2:
			# compute intersection of all optimal answer sets
			out_file_intersect = Path(outdir, f'results_intersect_{sep}.txt')
			flags_intersect = '--opt-mode=optN --quiet=1 --enum-mode cautious -W no-atom-undefined'
			show_file = Path(aspdir, 'show.pl')
			command_intersect = create_command(aspdir, [show_file, asp_input_file], out_file_intersect, sep, flags_intersect)
			subprocess.call(command_intersect, shell=True)

			# compute union of all optimal answer sets with the minimal cost from intersection computation,
			# saves computation of non-optimal answer sets.
			out_file_union = Path(outdir, f'results_union_{sep}.txt')
			cost = get_penalty_value_from_file(out_file_intersect)

			if cost != np.inf:
				flags_union = f'--opt-mode=enum,{cost} --quiet=1 --enum-mode brave -W no-atom-undefined'
				# Compute union
				command_union = create_command(aspdir, [show_file, asp_input_file], out_file_union, sep, flags_union)
				subprocess.call(command_union, shell=True)
			else:
				print('Did not compute union of all optimal answer sets. Computing the intersection of all optimal answer sets resultet in UNSAT.')

		elif mode == 3:
			get_graph_from_asp(asp_input_file, outdir, aspdir, sep, d)

		# save times
		t = timer() - start
		times.append({'sep': sep, 't': t})	
	
	# clean up
	files = [file for file in os.listdir(outdir) if re.match('^edge.*\.txt$', file) or re.match('^conf.*\.txt$', file)]
	[ os.remove(Path(outdir, file)) for file in files ]

	return times


def create_command(aspdir: Path, inputs: List[Path], out_file: Path, sep: str, flags: str) -> str:
	'''
	Creates string of a command that runs clingo.


	Returns:
	-------
	A string representing a command that calls clingo.

	Parameters:
	----------
	aspdir: Path
		Path to directory where ASP encodings live.
	inputs: List[Path]
		Paths of clingo input; CIRs files.
	out_file: Path
		Path to file where output should be saved.
	sep: str
		ASP encoding: s = σ-sep, d = d-sep(cyclic), a = d-sep(acyclic).
	flags: str
		A space separated string of clingo flags.
	# mode: int
	# 	1: add flags such that clingo finds one optimal answer set.
	# 	2: add flags such that clingo only computes and prints the optimization cost of one optimal answer set.
	# 	3: add flags such that clingo runs once per possible features with additional hard constraints.
	# intersect : bool
	# 	Only relevant for mode 2. If True compute intersect, else compute union.
	'''
	# flags
	# -W no-atom-undefined; is used to disable warning
	# --quiet=n; n can be 0 (means print all), 1 (means print last), and 2 (means do not print any models, costs, or individual call statistic)
	# --opt-mode=optN --quiet=1 computes all answer sets and only prints optimal answer sets
	# --enum-mode cautious computes intersection of all answer sets
	# --enum-mode brave computes intersection of all answer sets

	# choose encoding
	if sep == 's':
		encoding_file = "sigma_hej_cyclic.pl"
	elif sep == 'd':
		encoding_file = "hej_cyclic.pl"
	elif sep == 'a':
		encoding_file = "hej_acyclic.pl"

	inputs_str = ' '.join(str(p.absolute()) for p in inputs)

	command = (f'clingo {flags} ' 															     # flags
      		   f'{aspdir.absolute()}/{encoding_file} {aspdir.absolute()}/partial_comp_tree.pl '  # asp encoding
      		   f'{inputs_str}'                                                    				 # inputs
			   f' > {out_file.absolute()}' 														 # output
			  )

	return command


### Clingo output and input
def get_edge_from_asp(node0, node1, asp_input_file, filepath, aspdir,
					  edge_type='edge', sep='s'):
	# same logic as in original code
	# create temp files
	tmp_file_vs, tmp_file_pro = create_tmp_edge_files(node0, node1, filepath, edge_type)
	result_file_vs = Path(filepath, f'_sep_{sep}_results_vs.txt') 
	result_file_pro = Path(filepath, f'_sep_{sep}_results_pro.txt') 

	# create and call commands
	flags = '--quiet=2,1 -W no-atom-undefined'
	command_vs = create_command(aspdir, [asp_input_file, tmp_file_vs], result_file_vs, sep, flags)
	command_pro = create_command(aspdir, [asp_input_file, tmp_file_pro], result_file_pro, sep, flags)
	subprocess.call(command_vs, shell=True)
	subprocess.call(command_pro, shell=True)

	# parse output and remove tmps
	z, w = get_edge_from_files(result_file_vs, result_file_pro)
	os.remove(result_file_pro)
	os.remove(result_file_vs)
	return (z, w)


def get_graph_from_asp(asp_input_file: Path, filepath: Path, aspdir: Path, sep='s', d=5):
	# same logic as in original code
	edges = np.zeros(shape=(d, d))
	confs = np.zeros(shape=(d, d))
	edges_score = np.zeros(shape=(d, d))
	confs_score = np.zeros(shape=(d, d))

	for row in range(d):
		for col in range(d):  # edge, conf row <- col
			if col != row:
				(z_edge, w_edge) = get_edge_from_asp(col, row,
													 asp_input_file, filepath, aspdir,
													 edge_type='edge', sep=sep)
				edges[row, col] = z_edge
				edges_score[row, col] = w_edge

	edges_file_name = 'edges_pred_' + sep + '_sep.csv'
	np.savetxt(Path(filepath, edges_file_name).absolute(), edges, fmt='%1.1f', delimiter=",")
	edges_score_file_name = 'edges_score_' + sep + '_sep.csv'
	np.savetxt(Path(filepath, edges_score_file_name).absolute(), edges_score, fmt='%i', delimiter=",")
	
	for row in range(d):
		for col in range(d):
			if col > row:
				(z_conf, w_conf) = get_edge_from_asp(row, col,
													 asp_input_file, filepath, aspdir,
													 edge_type='conf', sep=sep)
				confs[row, col] = z_conf
				confs[col, row] = z_conf
				confs_score[row, col] = w_conf
				confs_score[col, row] = w_conf

	confs_file_name = 'confs_pred_' + sep + '_sep.csv'
	np.savetxt(Path(filepath, confs_file_name).absolute(), confs, fmt='%1.1f', delimiter=",")
	confs_score_file_name = 'confs_score_' + sep + '_sep.csv'
	np.savetxt(Path(filepath, confs_score_file_name).absolute(), confs_score, fmt='%i', delimiter=",")
	return (edges, confs, edges_score, confs_score)


def get_edge_from_files(file_vs, file_pro):
	# same logic as in original code
	w_vs = get_penalty_value_from_file(file_vs) #penalty against edge
	w_pro = get_penalty_value_from_file(file_pro) #penalty pro edge
	if w_vs == w_pro:
		confidence_pro_edge = 0
		z = 1/2
	else:
		confidence_pro_edge = w_vs - w_pro
		z = (1+np.sign(confidence_pro_edge))/2
	return (z,confidence_pro_edge) #absolute(confidence_pro_edge))


def create_tmp_edge_files(node0: int, node1: int, filepath: Path, edge_type: str='edge') -> None:
	# same logic as old code
	# create files for hard constraints
	tmp_file_vs = Path(filepath, f'{edge_type}_{node0}_{node1}_vs_.txt')
	tmp_file_pro = Path(filepath, f'{edge_type}_{node0}_{node1}_pro_.txt')

	with open(tmp_file_vs, 'w') as ff_vs:
		ff_vs.write(':-'+edge_type+'('+str(node0)+','+str(node1)+').')

	with open(tmp_file_pro, 'w') as ff_pro:
		ff_pro.write(edge_type+'('+str(node0)+','+str(node1)+').')

	return(tmp_file_vs, tmp_file_pro)


def get_penalty_value_from_file(path):
	# same logic as in original code
	go = True
	value = 0
	with open(path) as fl:
		for line in fl:
			if 'UNSATISFIABLE' in line:
				value = np.inf
				go = False
			elif go and 'SATISFIABLE' in line:
				value = 0
				go = False
			elif go and 'Optimization : ' in line:
				value = [int(s) for s in line.split() if s.isdigit()][0]
	return value


#### From conditional independence relations to asp constraints

def transform_tests_file_to_asp_file(df: pd.DataFrame, d: int, out_file: Path, use_CIT: bool = True, alpha: float = 0.01) -> None:
	'''
	Transforms conditional independece test file to asp file.


	Returns:
	-------
	The transformation is saved to a file. 

	Parameters:
	----------
	df : pd.DataFrame,
		Dataframe conditional independence relations (CIRs)
	d : int
		Number of observed variabels.
	out_file : Path
		Path to file where output should be saved. 
	use_CIT : bool
		If True, the expected format of the CIRs is: X, Y, do_targets, Z, p-val
		else: X, Y, do_targets, Z, sep where sep is True iff the nodes are seperated according to a markov property.
	'''
	asp_str = f"#const nrnodes = {d}. \n" ##const nrnodes = 4.

	for k in range(df.shape[0]):
		
		if use_CIT:
			weight, indep = transform_CIT(df['p-val'].iloc[k], al=alpha)
		else:
			weight, indep = 1, df['sep'].iloc[k]

		asp_str += strs_to_asp_expr(df['X'].iloc[k], df['Y'].iloc[k], df['Z'].iloc[k], df['do_targets'].iloc[k], d=d, w=weight, indep=indep)		 

	with open(out_file, 'a') as ff:
		ff.write(asp_str)


def transform_CIT(p_val, al=0.01, mul=1000, infty=1000):
	# same logic as old code; transforms p-value to weight and independence indicator
	# mul Multiplier for encoding p_values of independence tests to discrete values
	# infty Maximal  absolute value of weights (infty * mul)

	indep = (p_val>= al)
	if al != 0:
		odds = p_val/al
	else:
		odds = np.inf
	if odds <= 1e-316:
		score = -infty
	else:
		score = np.log(odds)
	if mul == np.inf:
		weight = infty
	else:
		weight =  np.minimum(mul*np.absolute(score),mul*infty)
	if np.isnan(weight):
		weight = 0
	weight = int(weight)
	
	return weight, indep


def list_to_number(lst):
	return sum({(2 ** f) for f in lst})


def str_to_list(s):
	return list({int(v) for v in s[1:-1].split(';') if v != ''})


def str_to_number(s):
	return list_to_number(str_to_list(s))


def strs_to_cmpl_number(X, Y, Z, do='[]', d=4):
	lst = str_to_list(X) + str_to_list(Y) + str_to_list(Z) + str_to_list(do)
	set1 = set(range(d)).difference(lst)
	return list_to_number(list(set1))


def strs_to_asp_expr(X, Y, Z, do='[]', d=4, w=0, indep=True):
	if indep:
		ind = 'in'
	else:
		ind = ''
	expr = ind + 'dep( %d , %d , %d , %d , %d , %d ). \n' % (
	str_to_list(X)[0], str_to_list(Y)[0], str_to_number(Z), str_to_number(do), strs_to_cmpl_number(X, Y, Z, '[]', d),
	int(w))
	return expr


#### Calculate conditinal independence relations

def run_all_CIT(data: List[Tuple[np.ndarray, List[int]]], d: int) -> pd.DataFrame:
	'''
	Runs all conditional independece test.


	Returns:
	-------
	A dataframe with results of conditional independence tests (CITs).
	Each row consits of X, Y, do_targets, Z, p-val and where p-val is the result of the CIT that
	tests if node X is independent of node Y given nodes Z with interventions on nodes in do_targets.   

	Parameters:
	----------
	data : List[np.ndarray, List[int]]
		An array of tuples. Each tuple (samples_i, interventions_i) consists of a n x d sample matrix,
		where n is the number of samples and d the number of observed variables, and an intervention vector
		that contains the indexes of all intervened variables for the given samples.
	d : int
		Number of observed variabels.
	'''
	partitions = all_xyZ_partitions_list(range(d))
	df = pd.DataFrame()
	i = 0

	for S, do_targets in data:
		
		S = normalranktransform(S)

		for id_x, id_y, id_z in partitions:
			df.loc[i,'X'] = '['+';'.join([str(s) for s in id_x])+']'
			df.loc[i,'Y'] = '['+';'.join([str(s) for s in id_y])+']'
			df.loc[i,'do_targets'] = '['+';'.join([str(s) for s in do_targets])+']'
			df.loc[i,'Z'] = '['+';'.join([str(s) for s in id_z])+']'
			df.loc[i,'p-val'] = pcor_p(S, id_x, id_y, id_z)
			i += 1

	return df

def powerlist(lst):
	pwlst = []
	d = len(lst)
	for i in range(2 ** d):
		subset = [x for j, x in enumerate(lst) if (i >> j) & 1]
		pwlst.append(subset)
	return pwlst


def all_xyZ_partitions_list(lst):
	'''
	Creates all partitions of a list into: two elements of the list x and y,
	and the rest of the the list Z. 

	Return:
	------
	A list of all unique partitions [[x], [y], Z]. 
	The list only inlucd one of [[x], [y], Z] and [[y], [x], Z],
	since they are considered the same partition.

	'''
	lst = set(lst)
	all_xyZ = []
	for i, x in enumerate(lst):
		for j, y in enumerate(lst):
			if j > i:
				ll = list(lst.difference({x, y}))
				for Z in powerlist(ll):
					all_xyZ.append(([x], [y], Z))
	return list(all_xyZ)


def normalranktransform(df):
	n = df.shape[0]
	d = df.shape[1]
	df = pd.DataFrame(df)
	# S = df.rank(method='first').as_matrix()/(n+1)
	S = df.rank(method='first').to_numpy()/(n+1)
	return scipy.stats.norm.ppf(S).reshape((n,d))


def pcor(df, id_x, id_y, id_z):
	ex = get_residuals(df, id_x, id_z)
	ey = get_residuals(df, id_y, id_z)
	# print(ey.shape, ex.shape)
	#reshape added; stats.pearsonr compile error
	ex = ex.reshape((ex.shape[0],))
	ey = ey.reshape((ey.shape[0],))
	pcr, dd = scipy.stats.pearsonr(ex, ey)
	return pcr


def p_value(pcr, n_samples, dim_cond_set):
	degfr = max(n_samples - dim_cond_set - 2, 1)
	pcr2 = pcr ** 2
	if pcr2 == 1:
		p_val = 0
	# elif degfr < 1:
	#    p_val = np.nan
	else:
		value = np.sqrt((degfr * pcr2) / (1. - pcr2))
		p_val = scipy.stats.t.sf(value, degfr) * 2
	return p_val


def pcor_p(df, id_x, id_y, id_z):
	pcr = pcor(df, id_x, id_y, id_z)
	n = df.shape[0]
	d = len(id_z)
	p_val = p_value(pcr, n, d)
	return p_val


def get_residuals(df, id_y, id_z):
	n = df.shape[0]
	Y = df[:, id_y]
	Y = Y - Y.mean(axis=0)  # .reshape(n,len(id_y))
	Y = Y / Y.std(axis=0)  # .reshape(n,len(id_y))
	Z = df[:, id_z]
	Z = Z - Z.mean(axis=0)  # .reshape(n,len(id_z))
	Z = Z / Z.std(axis=0)  # .reshape(n,len(id_z))
	if len(id_z) > 0:
		beta = np.linalg.lstsq(Z, Y, rcond=-1)[0]
		residual = Y - np.dot(Z, beta)
	else:
		residual = Y
	return residual
