# cancer_classification
Code accompanying "Uncertainty in Lung Cancer Stage for Outcome Estimation via Set-Valued Classification" [arXiv:2107.01251](https://arxiv.org/abs/2107.01251).

Simulation Study
1.	misc_funs_v3.R: contains underlying functions for executing simulation study 
2.	sim_server_set_1.R: simulate data, single iteration estimation, bootstrap estimation, all for scenario 1.
3.	sim_server_set_2.R: simulate data, single iteration estimation, bootstrap estimation, all for scenario 2.
4.	sim_server_set_3.R: simulate data, single iteration estimation, bootstrap estimation, all for scenario 3.
5.	sim_results_all.R: processes simulation results into table output

Data Analysis
1.	misc_funs_empirical_v3.R: contains underlying functions for executing data analysis 
2.	empirical_estimation.R: data analysis for naive standard practice, naive bootstrap, and weighted bootstrap
3.	empirical_results.R: converting results from empirical_estimation.R into tables and output

Combined
1.	results_figs.R: Converts tables / results from sim_results_all.R and empirical_estimation_results.R to manuscript figures (main text and appendix)
