using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,ProgressMeter
Random.seed!(123)
mainpath="/group/qtlchenggrp/tianjing/GWAS_wholepig/"
cd(mainpath)

include(mainpath*"GWAS4BayesNN.jl")

method="hmc_bayescpi_2"

mme=readdlm(mainpath*"mme.output_genotypes.txt",',')
map_file    = "pig_map_int.csv"  #map file for all chromosomes
window_size = "1 Mb"
threshold   = 0.01


weights_file=mainpath*"$method/results/MCMC_samples_neural_networks_bias_and_weights.txt"
marker_effects_file1=mainpath*"$method/results/MCMC_samples_marker_effects_geno_y1.txt"
marker_effects_file2=mainpath*"$method/results/MCMC_samples_marker_effects_geno_y2.txt"
ebv_file=mainpath*"$method/results/MCMC_samples_EBV_NonLinear.txt"
GWAS4BayesNN(mme,
	         map_filei,
	         marker_effects_file1,
	         marker_effects_file2,
	         weights_file,
	         ebv_file,
	         window_size = window_size,
	         threshold=thresholdi)