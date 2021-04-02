using Random, DelimitedFiles
Random.seed!(123)

cd("/home/tianjing/paper_simu/train_test_ID")

nObs=928
nRep =50


for i in 1:nRep
    ID      = randperm(nObs)
    
    trainID = ID[1:743]
    testID  = ID[744:end]
    
    open("rep$i.trainID.txt", "w") do io
        writedlm(io, trainID)
    end

    open("rep$i.testID.txt", "w") do io
        writedlm(io, testID)
    end
end
