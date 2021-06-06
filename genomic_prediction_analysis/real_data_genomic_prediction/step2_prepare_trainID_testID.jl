#this file is to generate trainID and testID for genomic prediction
using Random, DelimitedFiles, CSV, DataFrames
speciesALL=["maize","pig","soy","spruce","switchgrass"]
nRep = 20  #rep1,...,rep20

for species in speciesALL
    Random.seed!(123)
    cd("/group/qtlchenggrp/tianjing/paper_real/$species/train_test_ID")

    pheno_path="/home/tianjing/BNN_data/$species/$species"*"_pheno.csv"  #path for phenptype
    y=CSV.File(pheno_path) |> DataFrame   #read phenotype

    nObs = size(y,1)              #number of individuals
    nTrain = Int(floor(nObs*0.8)) #80% training, 20% testing

    for i in 1:nRep
        ID      = randperm(nObs)  #construct a random permutation of length nObs
        trainID = ID[1:nTrain]    #index for training set
        testID  = ID[nTrain+1:end]#index for testing set
        
        ### save train&test index.
        open("rep$i.trainID.txt", "w") do io
            writedlm(io, trainID)
        end

        open("rep$i.testID.txt", "w") do io
            writedlm(io, testID)
        end
    end
end






