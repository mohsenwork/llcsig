#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# arguments for calling function /llc/llc_bootstrap or /llc/llc_plain
# data_path     - path where pickled data is saved
# out_dir       - path to where output is saved
# mode          - specifies which version to run; 
#                 mode=1: plain version with finite samples. Runs LLC once. The argument mode should be 1
#                         The data in data_path must is pickled python dictionary with the intervention vector as keys 
#                         and sampls as values. The covariance matrix is calculated, then llc_plain is called.
#                 mode=2: bootstrapped version. Runs LLC with resampling. The argument mode should be 2.
#                         The data in data_path must is pickled python dictionary with the intervention vector as keys 
#                         and sampls as values. llc_bootstrap is called with the samples. 
#                 mode=3: plain version with oracle covariance. Runs LLC once. The argument mode should be 3
#                         The data in data_path must is pickled python dictionary with the intervention vector as keys 
#                         and oracle cov matrixes as values. llc_plain is called using the oracle covariance matrix.
# penalty       - 'none', 'L0', 'L1' or 'L2' - type of penalization/regularization used
#                  in the solving of the linear equation systems
# lambda        - regularization parameter
# n_bootstraps  - number times to run llc_plain (only relevant for mode 2)
# maxsetsize    - The maximum size of a set that is used to condition on when using faithfulness rules.
#                 -1 (or NA) for no restrictions.
# rules         - which of the faithfulness rules (1,2,3) should be used.
# null          - specifies whether to request identifiability information about the null-space. Use null=1 to request such information.



if (length(args) != 11)
  stop(paste0('Expecting 11 arguments, got ', length(args)), call.=FALSE)

# parse input parameters
data_path = args[1]
out_dir = args[2]
mode = as.numeric(args[3])
penalty = args[4]
lambda = as.numeric(args[5])
n_bootstraps = as.numeric(args[6])
maxsetsize = as.numeric(args[7])
tmp = args[8]
rules = as.numeric(unlist(strsplit(tmp, '')))
if ( any(rules == 0) ) {
  rules = c()
}
null = as.logical(as.numeric(args[9]))
alpha = as.numeric(args[10])
verbose = as.logical(as.numeric(args[11]))

# load code
library(reticulate)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '\\.[RrSsQq]$')) {
    if(trace) cat(nm,':')
    source(file.path(path, nm), ...)
    if(trace) cat('\n')
  }
}

if (verbose)
  cat(paste0('INFO - ', Sys.time(), ': Loading the R-code...\n'))
sourceDir('./simulations/resources/llc/2012/model/',trace=FALSE)
sourceDir('./simulations/resources/llc/2012/llc/',trace=FALSE)
sourceDir('./simulations/resources/llc/2012/common/',trace=FALSE)

if (verbose)
  cat(paste0('INFO - ', Sys.time(), ': Loading data and running LLC...\n'))

# parse data
df = py_load_object(data_path, pickle = 'pickle')
Data = list()
i = 1

for (name in names(df)) {

    D = list()
    sub = substr(name, 2, nchar(name)-1)
    interventions = as.numeric(unlist(strsplit(sub, ', ')))
    D$e = interventions
    # Number of variables
    n = length(interventions)

    if (mode == 1){
      # input is data for run plain
      data = t(df[[name]])
      D$Cx = cov(data)
      # Weight for each experiment is number of samples
      D$N = dim(df[[name]])[2]

    } else if (mode == 2){
      # input is data for run plain
      D$data = t(df[[name]])
      D$N = dim(df[[name]])[2]

    } else if (mode == 3){
      # input is oracle cov for run plain in inf sample limit
      data = t(df[[name]])
      D$Cx = data
      # Weight for each experiment is 1
      D$N = Inf
    }
    Data[[i]] = D
    i = i+1
}

# run llc
if ( length(rules) > 0) {
  if (verbose)
    cat(paste0('INFO - ', Sys.time(), ': Running LLC with faithfulness rules ..\n'))
} else {
  if (verbose)
    cat(paste0('INFO - ', Sys.time(), ': Running LLC without faithfulness rules ..\n'))  
}

if (mode == 1) {
  if (verbose)
    cat(paste0('INFO - ', Sys.time(), ': Running llc_plain with data..\n'))
  out = llc_plain(DATA=Data, maxsetsize=maxsetsize, rules=rules, penalty=penalty, lambda=lambda, null=null, alpha=alpha)
}else if (mode == 2){
  if (verbose)
    cat(paste0('INFO - ', Sys.time(), ': Running llc_bootstrap with data..\n'))
  out = llc_bootstrap(DATA=Data, n=n, maxsetsize=maxsetsize, rules=rules, penalty=penalty, lambda=lambda, alpha=alpha, bootstraps=n_bootstraps)
}else if (mode == 3){
  if (verbose)
    cat(paste0('INFO - ', Sys.time(), ': Running llc_plain with oracle input..\n'))
  out = llc_plain(DATA=Data, maxsetsize=maxsetsize, rules=rules, penalty=penalty, lambda=lambda, null=null, alpha=alpha)
}

# save output
py_save_object(out, paste0(out_dir, '/llc_out.pkl'), pickle = 'pickle')
if (verbose)
  cat(paste0('INFO - ', Sys.time(), ': Done.\n'))