# this file is the Julia scripts to calculate prediction accuracy
ENV["GKSwstype"]="100"  #to enable saving figure in server
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles,Plots

species     = ARGS[1]
method      = ARGS[2]
nIter       = 20_000
nrep        = 20  #rep1,...,rep20

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

trait_symbol = Symbol(traitname)  #e.g, :YLD

mainpath="/group/qtlchenggrp/tianjing/paper_real/$species/"
cd(mainpath*"$method/")

############ READ DATA ##########
phenofile="/home/tianjing/BNN_data/$species/$species"*"_pheno.csv" #path for phenotype
phenotypes = DataFrame!(CSV.File(phenofile));                      #read phenptype
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);                      #set individuals ID as string


for rep in 1:nrep  #rep1,..,rep20
    #read individuals ID for testing
    id_path = mainpath*"train_test_ID/"
    testID = vec(Int64.(readdlm(id_path*"rep$rep.testID.txt")))
    
    #read MCMC samples for EBV
    if occursin("hmc",method)  #NN-Bayes
        out_ebv = CSV.read("rep$rep/results/EBV_NonLinear.txt",delim = ',',header=true,DataFrame)
        MCMC_ebv = CSV.read("rep$rep/results/MCMC_samples_EBV_NonLinear.txt",delim = ',',header=true,DataFrame)
    elseif !occursin("hmc",method) #linear model
        out_ebv = CSV.read("rep$rep/results/EBV_$traitname.txt",delim = ',',header=true,DataFrame)
        MCMC_ebv = CSV.read("rep$rep/results/MCMC_samples_EBV_$traitname.txt",delim = ',',header=true,DataFrame)
    end
    
    out_ebv[!,:ID]=string.(out_ebv[!,:ID]); #get individual ID
    MCMC_ebv = Matrix(MCMC_ebv)'            #a matrix of (n,1000), where n is number of individuals, 1000 is MCMC samples
    #For example, for chain of length 20000, MCMC chain was saved for every 20 iterstions, to form 1000 MCMC samples.
    #the number of samples can be changed via runMCMC(output_samples_frequency=)

    #calculate mean of 1000 MCMC samples
    cum_ebv = cumsum(MCMC_ebv,dims=2)
    for i in 1:1000  
        cum_ebv[:,i]=cum_ebv[:,i]/i
    end
    cum_ebv=DataFrame(cum_ebv)
    insertcols!(cum_ebv, 1, :ID => out_ebv[!,:ID])

    #combine true phenotype and ebv
    results = innerjoin(cum_ebv, phenotypes, on = :ID)

    #calculate prediction accuracy for each MCMC samples
    accuracy=zeros(1000)
    for i in 1:1000
        accuracy[i] = cor(results[testID,Symbol("x$i")],results[testID,trait_symbol])
    end

    #same results
    open("accuracy.$method.rep$rep.txt", "w") do io
        writedlm(io, accuracy)
    end

    #visualization
    myfig=plot(collect(nIter/1000:nIter/1000:nIter),accuracy)
    savefig(myfig,"accuracy.$method.rep$rep.png")
end

