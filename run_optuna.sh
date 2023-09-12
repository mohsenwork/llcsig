#!/bin/bash
make run_optuna setting=1 name=optuna_llc_nf_l1_roc conf=optuna_llc_nf_l1.yml n_trials=300
make run_optuna setting=1 name=optuna_llc_nf_l2_roc conf=optuna_llc_nf_l2.yml n_trials=300
make run_optuna setting=2 name=optuna_llc_nf_l1_acc conf=optuna_llc_nf_l1.yml n_trials=300
make run_optuna setting=2 name=optuna_llc_nf_l2_acc conf=optuna_llc_nf_l2.yml n_trials=300
make run_optuna setting=3 name=optuna_llc_f_l1_roc conf=optuna_llc_f_l1.yml n_trials=300
make run_optuna setting=3 name=optuna_llc_f_l2_roc conf=optuna_llc_f_l2.yml n_trials=300
make run_optuna setting=4 name=optuna_llc_f_l1_acc conf=optuna_llc_f_l1.yml n_trials=300
make run_optuna setting=4 name=optuna_llc_f_l2_acc conf=optuna_llc_f_l2.yml n_trials=300
make run_optuna setting=5 name=optuna_asp_d_roc conf=optuna_asp_d.yml n_trials=300
make run_optuna setting=6 name=optuna_asp_d_acc conf=optuna_asp_d.yml n_trials=300
make run_optuna setting=5 name=optuna_asp_s_roc conf=optuna_asp_s.yml n_trials=300
make run_optuna setting=6 name=optuna_asp_s_acc conf=optuna_asp_s.yml n_trials=300
make run_optuna setting=7 name=optuna_asp_d_acc_narrow conf=optuna_asp_d.yml n_trials=300
make run_optuna setting=7 name=optuna_asp_s_acc_narrow conf=optuna_asp_s.yml n_trials=300
