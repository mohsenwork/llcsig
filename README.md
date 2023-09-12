# Cyclic Causality

## Installation

To run this project, the following packages are required:

- Python
- R
- clingo

### Python Packages

The Python packages required for this project are listed below, along with the tested versions:

- pyyaml==6.0
- numpy==1.22.3
- pandas==1.4.4
- scikit-learn==1.1.1
- seaborn==0.12.2
- matplotlib==3.6.3
- optuna==3.1.0
- networkx==2.8.4
- scipy==1.10.0
- pytest==7.1.2
- pytest-cov==4.0.0

To install Python and the required packages, please follow the instructions below:

1. Install Conda (for Linux systems):
   - Download the Anaconda installer script specifically for Linux:

     ```shell
     wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
     ```

   - Run the installer script:

     ```shell
     bash Anaconda3-2022.10-Linux-x86_64.sh
     ```

2. Create a Conda environment and install Python packages:

   - Add conda-forge channel to Conda:

      ```shell
      conda config --add channels conda-forge 
      ```
   - Create a Conda environment named `cycau-env`:

     ```shell
     conda create -n cycau-env python=3.10.8
     ```

   - Activate the Conda environment:

     ```shell
     conda activate cycau-env
     ```

   - Install the required Python packages using Conda and the specified versions:

     ```shell
     conda install -n cycau-env -c plotly plotly==5.13.1 pyyaml==6.0 numpy==1.22.3 pandas==1.4.4 seaborn==0.12.2 matplotlib==3.6.3 optuna==3.1.0 networkx==2.8.4 scipy==1.10.0 pytest==7.1.2 pytest-cov==4.0.0 scikit-learn==1.1.1
     ```
     
### R Packages

To install R and the required R package, please follow the instructions below:

1. Install R:
   - On Ubuntu, run the following command to install R:

     ```shell
     sudo apt-get update && apt-get install -y r-base
     ```

2. Install the R package:
   - Open an R console and install the required R package using the following command:

     ```R
     install.packages(c('reticulate'), repos='https://cran.rstudio.com/')
     ```

### Clingo

Clingo can be installed using Conda with the following command:

```shell
conda install -n cycau-env -c potassco clingo
```

### Jupyter Notebook 
Jupyter notebook can be installed with the following command:

```shell
conda install -c conda-forge notebook==6.4.12
```

## Running the Code

### Reproducing Figures 1-8 in the Paper

To reproduce Figures 1-8 in the paper, execute the following command:

```shell
make run_multiprocessing conf=hyperparameter_eval.yml n_models=500
```

- The above command will run simulations based on the configuration specified in the `hyperparameter_eval.yml` file.
- The simulations will be executed using multiprocessing to speed up the process.
- The resulting data will be saved in the "data" subdirectory.

After running the above command, the data required to reproduce the figures will be saved. To generate the figures, you can use the provided Jupyter notebook named `hyperparameter_eval.ipynb` in the `plots`  directory. Open the notebook and follow the instructions within to load the saved data and reproduce Figures 1-8.


### Running the Optuna Studies

The Optuna studies in this project determine the hyperparameters as discussed in the paper. These studies can take a long time to run (!). 

To run all the studies listed below, you can use the script `run_optuna.sh`.
Make it executable using the command: 
```shell
chmod +x run_optuna.sh
```

Then, you can run the script:
```shell
./run_optuna.sh
```

#### LLC Studies

The following studies determine the hyperparameters for LLC-NF and LLC-F algorithms as discussed in the paper. To visualize the studies, use the notebook `llc_optuna.ipynb` located in the `plots` directory. However, before visualizing the studies, we need to generate the data, which will be saved in the `plots/optuna` directory. To run the studies for LLC algorithms and generate the necessary data, follow the instructions below:

#### LLC-NF Algorithm Studies

1. Study using AUC ROC as metric and l1 penalty:
   ```shell
   make run_optuna setting=1 name=optuna_llc_nf_l1_roc conf=optuna_llc_nf_l1.yml n_trials=300
   ```

2. Study using AUC ROC as metric and l2 penalty:
   ```shell
   make run_optuna setting=1 name=optuna_llc_nf_l2_roc conf=optuna_llc_nf_l2.yml n_trials=300
   ```

3. Study using accuracy as metric and l1 penalty:
   ```shell
   make run_optuna setting=2 name=optuna_llc_nf_l1_acc conf=optuna_llc_nf_l1.yml n_trials=300
   ```

4. Study using accuracy as metric and l2 penalty:
   ```shell
   make run_optuna setting=2 name=optuna_llc_nf_l2_acc conf=optuna_llc_nf_l2.yml n_trials=300
   ```

#### LLC-F Algorithm Studies

1. Study using AUC ROC as metric and l1 penalty:
   ```shell
   make run_optuna setting=3 name=optuna_llc_f_l1_roc conf=optuna_llc_f_l1.yml n_trials=300
   ```

2. Study using AUC ROC as metric and l2 penalty:
   ```shell
   make run_optuna setting=3 name=optuna_llc_f_l2_roc conf=optuna_llc_f_l2.yml n_trials=300
   ```

3. Study using accuracy as metric and l1 penalty:
   ```shell
   make run_optuna setting=4 name=optuna_llc_f_l1_acc conf=optuna_llc_f_l1.yml n_trials=300
   ```

4. Study using accuracy as metric and l2 penalty:
   ```shell
   make run_optuna setting=4 name=optuna_llc_f_l2_acc conf=optuna_llc_f_l2.yml n_trials=300
   ```

### ASP Studies

The Optuna studies in this project also include ASP algorithm studies. These studies determine the hyperparameters for the ASP-d and ASP-s algorithms. 
To visualize the studies, use the notebook `asp_optuna.ipynb` located in the `plots` directory. However, before visualizing the studies, we need to generate the data, which will be saved in the `plots/optuna` directory. To run the studies for LLC algorithms and generate the necessary data, follow the instructions below:

To run the ASP studies, follow the instructions below:

#### ASP-d Algorithm Studies

1. Study using AUC ROC as metric:
   ```shell
   make run_optuna setting=5 name=optuna_asp_d_roc conf=optuna_asp_d.yml n_trials=300
   ```

2. Study using accuracy as metric:
   ```shell
   make run_optuna setting=6 name=optuna_asp_d_acc conf=optuna_asp_d.yml n_trials=300
   ```

#### ASP-s Algorithm Studies

1. Study using AUC ROC as metric:
   ```shell
   make run_optuna setting=5 name=optuna_asp_s_roc conf=optuna_asp_s.yml n_trials=300
   ```

2. Study using accuracy as metric:
   ```shell
   make run_optuna setting=6 name=optuna_asp_s_acc conf=optuna_asp_s.yml n_trials=300
   ```

#### Narrow Search Space Studies

To perform ASP studies with a narrower search space,

 follow the instructions below:

#### ASP-d Algorithm Narrow Search Space Study

1. Study using accuracy as metric with a narrower search space:
   ```shell
   make run_optuna setting=7 name=optuna_asp_d_acc_narrow conf=optuna_asp_d.yml n_trials=300
   ```

#### ASP-s Algorithm Narrow Search Space Study

1. Study using accuracy as metric with a narrower search space:
   ```shell
   make run_optuna setting=7 name=optuna_asp_s_acc_narrow conf=optuna_asp_s.yml n_trials=300
   ```

For any issues or questions, please contact the project maintainers.

## Project structure
 ```
├── config                                                         # Configuration files containt values for arguments when calling functions
│   └──  ...
├── Makefile
├── plots                                                          # Directory for notebooks
│   └──  ...
├── pytest.ini
├── README.md
├── setup.cfg
├── simulations
│   ├── app.py                                                     # Entry point for running ASP and LLC
│   ├── common                                                     # Utility functions
│   │   └── ...
│   ├── data_generation                                            # Code for creating data from models
│   │   └──  ...
│   ├── __init__.py
│   ├── __main__.py                                                # Main
│   ├── multiprocess_simulations.py                                # Code for multiprocess running
│   ├── objective_1.py                                             # Optuna objective files - start
│   ├── objective_2.py
│   ├── objective_3.py
│   ├── objective_4.py
│   ├── objective_5.py
│   ├── objective_6.py
│   ├── objective_7.py
│   ├── objective.py
│   ├── optuna_custome.py                                         # Optuna objective files - end
│   └── resources
│       ├── llc                                                    # LLC code
│       │   └── 2012
│       │       ├── common
│       │       │   └── ...
│       │       ├── llc
│       │       │   └── ...
│       │       ├── main.R                                         # Entry point for running LLC code
│       │       └── model
│       │           └── ...
│       └── sigmasep                                               # ASP code
│           ├── ASP
│           │   └── ...
│           └── main.py                                            # Entry point for running ASP code
└── tests                                                          # Testing directory
    └── ...
```
## Tested Environment

The code for this project has been tested on an AWS t3.large instance with the following specifications:

- Instance Type: AWS t3.large
- RAM: 8GB
- Storage: 24GB
- Operating System: Ubuntu 22.04 LTS

Please note that while this is the tested environment, the project also works on other similar configurations.

If you encounter the error message "mkl-service + Intel(R) MKL: MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library" on AWS instances, you can resolve it by exporting MKL_THREADING_LAYER=GNU. This command ensures compatibility between the mkl-service library and libgomp.so.1.