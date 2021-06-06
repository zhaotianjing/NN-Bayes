# this file is the Julia scripts to calculate prediction accuracy

ENV["GKSwstype"]="100" #to enable saving figure in server
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles,Plots

repy              = parse(Int, ARGS[1])
method            = ARGS[2]
chainLength       = 10_000
nrep              = 50  #rep1,...,rep50
phenoname="pig.y.epistasis.rep$repy.txt"

mainpath="/home/tianjing/paper_simu/"
phenofile  = mainpath*"y$repy/pig.y.epistasis.rep$repy.txt"

cd(mainpath*"y$repy/$method/")

############ READ DATA ##########
#read phenotype
phenotypes = DataFrame!(CSV.File(phenofile));
phenotypes[!,:ID]=string.(phenotypes[!,:ID]);

for rep in 1:nrep  #rep1,..,rep50
    ##read individuals ID for testing
    id_path = mainpath*"train_test_ID/"
    testID = vec(Int64.(readdlm(id_path*"rep$rep.testID.txt")))
    
    #read MCMC samples for EBV
    if occursin("hmc",method)  #NN-Bayes
        out_ebv_y = CSV.read("rep$rep/results/EBV_NonLinear.txt",delim = ',',header=true,DataFrame)
        MCMC_ebv_y = CSV.read("rep$rep/results/MCMC_samples_EBV_NonLinear.txt",delim = ',',header=true,DataFrame)
    elseif !occursin("hmc",method)  #linear models
        out_ebv_y = CSV.read("rep$rep/results/EBV_y.txt",delim = ',',header=true,DataFrame)
        MCMC_ebv_y = CSV.read("rep$rep/results/MCMC_samples_EBV_y.txt",delim = ',',header=true,DataFrame)
    end
        
    out_ebv_y[!,:ID]=string.(out_ebv_y[!,:ID]);  #get individual ID
    MCMC_ebv_y = Matrix(MCMC_ebv_y)'  #a matrix of (n,1000), where n is number of individuals, 1000 is MCMC samples
    #For example, for chain of length 20000, MCMC chain was saved for every 20 iterstions, to form 1000 MCMC samples.
    #the number of samples can be changed via runMCMC(output_samples_frequency=)

    #calculate mean of 1000 MCMC samples
    cum_ebv_y = cumsum(MCMC_ebv_y,dims=2)
    for i in 1:1000  
        cum_ebv_y[:,i]=cum_ebv_y[:,i]/i
    end
    cum_ebv_y=DataFrame(cum_ebv_y) 
    insertcols!(cum_ebv_y, 1, :ID => out_ebv_y[!,:ID])

    #combine true phenotype and ebv
    results_y = innerjoin(cum_ebv_y, phenotypes, on = :ID)

    #calculate prediction accuracy for each MCMC samples
    accuracy_y=zeros(1000)
    for i in 1:1000
        accuracy_y[i] = cor(results_y[testID,Symbol("x$i")],results_y[testID,:y])
    end

    #same results
    open("accuracy.y$repy.rep$rep.$method.txt", "w") do io
        writedlm(io, accuracy_y)
    end

    #visualization
    myfig=plot(collect(chainLength/1000:chainLength/1000:chainLength),accuracy_y,xlabel="iteration",ylabel="accuracy",title="$method.y$repy")
    savefig(myfig,"accuracy.y$repy.rep$rep.$method.png")
end

