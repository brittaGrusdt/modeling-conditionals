# modeling-conditionals

This repository contains all code and simulation data for our submission to S+P titled **Probabilistic Modeling of Rational Communication with Conditionals**.

All configuration data (parameters etc.) is stored in "config.yml".

To run the model with the default context prior, run "R/default-model/run-default-model.R". Specify on top of this file whether you want to run a sweep over parameter space (*run_sweep*), the name of the *target*-folder, a random or specific seed to run wppl model (*seed_wppl*) and whether the results shall be stored to a subfolder with the seed as name (*save_seed_subdir*). The data will be stored in this folder as subdirectory of "data/default-model".
With the configuration as it is set, the script will run once to get predictions of the pragmatic listener (given utterance 'If A, C'), once to get speaker predictions for 10,000 states sampled from the prior, once for the literal-speaker condition, i.e., for 10,000 states sampled from the prior, but for which the utterance 'If A, C' is assertable and once it returns 10,000 states sampled from the prior.

To run the model with one of the contexts of the Douven's cases, run "R/douven-examples/run-and-plot-douven-examples.R". 
This script runs each of the three cases presented in the paper (garden party, skiing and sundowners). 
The data will be stored in the your specified target folder as subdirectory of "data/douven-examples".


