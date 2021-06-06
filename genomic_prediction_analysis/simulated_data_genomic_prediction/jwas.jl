# this file is the Julia scripts to run JWAS for simulated datasets
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles
Random.seed!(123)

repy              = parse(Int, ARGS[1])  #y1,....,y20
rep               = parse(Int, ARGS[2])  #rep1,...,rep50
method            = ARGS[3]
method_cpy        = method
chainLength       = 10000

mainpath  ="/home/tianjing/paper_simu/"
genofile  = mainpath*"pig.x.n928.p5024.txt"     #genetype of chromosome1 of pig
phenofile = mainpath*"y$repy/pig.y.epistasis.rep$repy.txt"

cd(mainpath*"y$repy/$method/rep$rep/")

#obtain number of hidden nodes
if occursin("hmc",method)
	num_latent_traits = parse(Int,split(method,"_")[3])  #e.g., 2
	method=split(method,"_")[2]  # e.g., bayescpi
end

#set JWAS parameters
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
# read phenotype
phenotypes = DataFrame!(CSV.File(phenofile));
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);

# read trainID
id_path = mainpath*"train_test_ID/"
trainIDs = vec(Int64.(readdlm(id_path*"rep$rep.trainID.txt")))

############ JWAS ##########
geno = get_genotypes(genofile,method=jwas_method,estimatePi=jwas_estimatePi);

model_equations = "y = intercept + geno"; #set model equation

if occursin("hmc",method_cpy)       # NN-Bayes
	model = build_model(model_equations,num_latent_traits=num_latent_traits,nonlinear_function="Neural Network",activation_function="tanh");
	out   = runMCMC(model,phenotypes[trainIDs,:],mega_trait=true,chain_length=chainLength,double_precision=double_precision);
elseif !occursin("hmc",method_cpy)  # linear models
	model = build_model(model_equations);
	out   = runMCMC(model,phenotypes[trainIDs,:],chain_length=chainLength,double_precision=double_precision);
end







