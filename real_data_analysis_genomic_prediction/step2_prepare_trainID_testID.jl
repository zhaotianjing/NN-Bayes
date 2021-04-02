#this file is to generate trainID and testID for genomic prediction
using Random, DelimitedFiles, CSV, DataFrames
speciesALL=["maize","pig","soy","spruce","switchgrass"]
nRep = 20

for species in speciesALL
    Random.seed!(123)
    cd("/group/qtlchenggrp/tianjing/paper_real/$species/train_test_ID")

    pheno_path="/home/tianjing/BNN_data/$species/$species"*"_pheno.csv"
    y=CSV.File(pheno_path) |> DataFrame

    nObs = size(y,1)
    nTrain = Int(floor(nObs*0.8))

    for i in 1:nRep
        ID      = randperm(nObs)
        trainID = ID[1:nTrain]
        testID  = ID[nTrain+1:end]
        
        open("rep$i.trainID.txt", "w") do io
            writedlm(io, trainID)
        end

        open("rep$i.testID.txt", "w") do io
            writedlm(io, testID)
        end
    end
end






