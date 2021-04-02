using DataFrames,CSV,Random,JWAS,DelimitedFiles,Statistics,ProgressMeter
include("/group/qtlchenggrp/tianjing/GWAS/GWAS4BayesNN.jl")
Random.seed!(123)

repy   = parse(Int, ARGS[1])
method = ARGS[2]

main_path="/group/qtlchenggrp/tianjing/GWAS/"
cd(main_path*"y$repy/$method")

mme      = readdlm(main_path*"mme.output_genotypes.txt",',')
map_file = main_path*"pig_map_chr1.csv"

window_size = "1 Mb"
threshold = 0.01
@show repy,method,threshold


if occursin("hmc",method)   #hmc
	weights_file="results/MCMC_samples_neural_networks_bias_and_weights.txt"
	marker_effects_file1="results/MCMC_samples_marker_effects_geno_y1.txt"
	marker_effects_file2="results/MCMC_samples_marker_effects_geno_y2.txt"
	ebv_file="results/MCMC_samples_EBV_NonLinear.txt"
	GWAS4BayesNN(mme,
                 map_file,
             	 marker_effects_file1,
             	 marker_effects_file2,
             	 weights_file,
             	 ebv_file,
             	 window_size = window_size,
             	 threshold=threshold)
elseif !occursin("hmc",method)  #linear model
	marker_effects_file="results/MCMC_samples_marker_effects_geno_y.txt"
	GWAS(mme,map_file,marker_effects_file,threshold=threshold,window_size=window_size)
end
