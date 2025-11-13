# Dynamics_of_Trachoma_SI
Data and code for "Dynamics of trachoma infection in West Africa revealed by a hidden state model"

Gambia_Data.csv: Contains individual level metadata and the longitudinal test results.

Trachoma_Main.R: Source this file to undertake the analysis.

Trachoma_Main_Continue.R: Loads a previous MCMC trace and continues. Useful for longer chains.

Trachoma_MCMC_rif.R: Gibbs updates for observation parameters and initial state parameters.

Trachoma_cpp_rif.cpp: Includes to IFFBS update, and functions for evaluating the likelihoods of Metropolis Hastings updates.
