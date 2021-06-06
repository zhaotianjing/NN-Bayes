#this file is the script to run JWAS.

using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles
Random.seed!(123)

repy              = parse(Int, ARGS[1])  #y1,...,y20
method            = ARGS[2]
method_cpy        = method
chainLength       = 500_000
burnin            = 200_000
output_samples_frequency = 50

mainpath="/group/qtlchenggrp/tianjing/GWAS/"
cd(mainpath*"y$repy/$method/")

traitname="y"
phenoname="pig.y.epistasis.rep$repy.txt"
genofile = mainpath*"pig.x.n928.p5024.txt"   #chromosome1
phenofile  = mainpath*"y$repy/"*phenoname

#obtain number of hidden nodes
if occursin("hmc",method) #NN-Bayes
	num_latent_traits = parse(Int,split(method,"_")[3]) #e.g., 2
	method=split(method,"_")[2]   #e.g., bayescpi
end

#set JWAS parameters for each method
if method=="rrblup"
	jwas_method="RR-BLUP"
	jwas_estimatePi=false
	double_precision=false
elseif method=="bayesa"
	jwas_method="BayesA"
	jwas_estimatePi=false
	double_precision=true
elseif method=="bayesb"
	jwas_method="BayesB"
	jwas_estimatePi=true
	double_precision=true
elseif method=="bayescpi"
	jwas_method="BayesC"
	jwas_estimatePi=true
	double_precision=false
elseif method=="bayesl"
	jwas_method="BayesL"
	jwas_estimatePi=false
	double_precision=false
else 
	println("error method")
end

############ READ DATA ##########
#read phenotype
phenotypes = DataFrame!(CSV.File(phenofile));
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);

geno = get_genotypes(genofile,method=jwas_method,estimatePi=jwas_estimatePi);
model_equations = "$traitname = intercept + geno";


############ RUN JWAS ##########
if occursin("hmc",method_cpy) #NN-Bayes
	model = build_model(model_equations,num_latent_traits=num_latent_traits,nonlinear_function="Neural Network",activation_function="tanh");
	out   = runMCMC(model,phenotypes,mega_trait=true,chain_length=chainLength,burnin=burnin,output_samples_frequency=output_samples_frequency,double_precision=double_precision);
elseif !occursin("hmc",method_cpy) #linear model
	model = build_model(model_equations);
	out   = runMCMC(model,phenotypes,chain_length=chainLength,burnin=burnin,output_samples_frequency=output_samples_frequency,double_precision=double_precision);
end

