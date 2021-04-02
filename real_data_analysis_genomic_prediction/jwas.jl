using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles
Random.seed!(123)


species           = ARGS[1]
method            = ARGS[2]
rep               = parse(Int, ARGS[3])  #rep1,...,rep20
chainLength       = 20_000 
method_cpy        = method

if species=="soy"
	traitname="YLD"
elseif species=="maize"
	traitname="FT"
elseif species=="spruce"
	traitname="DE"
elseif species=="switchgrass"
	traitname="AN"
elseif species=="pig"
	traitname="y"
else
	println("error")
end

mainpath="/group/qtlchenggrp/tianjing/paper_real/$species/"
cd(mainpath*"$method/rep$rep/")

############ READ DATA ##########
#read phenotype
phenofile="/home/tianjing/BNN_data/$species/$species"*"_pheno.csv"
phenotypes = DataFrame!(CSV.File(phenofile));
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);
#read genotype
genofile = "/home/tianjing/BNN_data/$species/$species"*"_geno_clean.txt"
# read trainID for one specific CV
id_path = mainpath*"train_test_ID/"
trainIDs = vec(Int64.(readdlm(id_path*"rep$rep.trainID.txt")))


if occursin("hmc",method)
	num_latent_traits = parse(Int,split(method,"_")[3])  #e.g., 2
	method=split(method,"_")[2]  # e.g., bayescpi
end

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


geno = get_genotypes(genofile,method=jwas_method,estimatePi=jwas_estimatePi);
model_equations = "$traitname = intercept + geno";

if occursin("hmc",method_cpy)   #NN-Bayes
	model = build_model(model_equations,num_latent_traits=num_latent_traits,nonlinear_function="Neural Network");
	out   = runMCMC(model,phenotypes[trainIDs,:],mega_trait=true,chain_length=chainLength,double_precision=double_precision);
elseif !occursin("hmc",method_cpy)  #linear models
	model = build_model(model_equations);
	out   = runMCMC(model,phenotypes[trainIDs,:],chain_length=chainLength,double_precision=double_precision);
end






