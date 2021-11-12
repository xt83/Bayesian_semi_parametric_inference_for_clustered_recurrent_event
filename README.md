# Bayesian_semi_parametric_inference_for_clustered_recurrent_event_with_a_zero_inflation_and_a_terminal_event
The current folder contains files for implementing the BMZ-DP approach introduced in "Bayesian semi-parametric inference for clustered recurrent events with zero-inflation and a terminal event"

we develop a Bayesian semi-parametric model to jointly characterize the zero-inflated recurrent event process and the terminal event process. We use a point mass mixture of non-homogeneous Poisson processes to describe the recurrent intensity and introduce shared random effects from different sources to bridge the non-terminal and terminal event processes. To enable robustness, we consider nonparametric Dirichlet processes to model the residual of the accelerated failure time model for the survival process as well as the cluster-specific frailty distribution, and we develop a Markov Chain Monte Carlo algorithm to enable efficient posterior inference.


List of Files:
1) combn_sim_s1_bmz_dp.R: Simulation study for proposed BMZ-DP model with piecewise constant function for baseline hazard of recurrent process
2) combn_sim_s2_bmz_dp.R: Simulation study for proposed BMZ-DP model with Weibull distribution for baseline hazard of recurrent process
3) real_data_analysis_wp.R: Analysis of real STRIDE data based on BMZ-DP with piecewise constant function for baseline hazard of recurrent process

Notes: 
1) The sample size and cluster number is based on n and num_hos
2) You will need to change path/directory names before running the example program. 
3) The codes for simulation study contains the data generating process. 
