include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)

testSim = cellSimulation(testGS, testIC, 1000)

#examine information content different parts of the simulation provide about other parts
#across different simulations, of different sizes
#use mutual information and transfer entropy
using TransferEntropy
using Statistics
#create data set with many cells from the same gene set
mcds = makeMultiCellDataSet(testGS, testIC, 500, 2000)

mcdsRNA = Array{Int}(undef, 10, 500, 2000)
for i=1:2000
    mcdsRNA[:,:,i] = mcds[i].SplicedRNA
end

#mutual information between gene pairs in each cell
est =VisitationFrequency(RectangularBinning(200))
miRes = Array{Float64}(undef, 10, 10, 2000)
for i=1:2000
    for j=1:10
        for k=1:10
            miRes[j,k,i] = mutualinfo(mcdsRNA[j,:,i], mcdsRNA[k,:,i], est)
        end
    end
end

#sum matrices to get total MI for each gene pair
miResSum = sum(miRes, dims=3)[:,:,1]
#zero the diagonal
for i=1:10
    for j=1:10
        if i==j
            miResSum[i,j] = 0
        end
    end
end

geneTFAMat = getTFGeneAMat(testGS)
rankVec = Vector{Vector{Int}}(undef, 10)
for i=1:10
    causalGenes = findall(x->x != 0, geneTFAMat[:,i])
    ranks = findall(in(causalGenes), reverse(sortperm(miResSum[:,i])))
    rankVec[i] = ranks
end

#conditional MI
cmiRes = Array{Float64}(undef, 10, 10, 2000)
for i=1:2000
    for j=1:10
        for k=1:10
            cmiRes[j,k,i] = conditional_mutualinfo(mcdsRNA[j,:,i], mcdsRNA[k,:,i], est)
        end
    end
end

#sum matrices to get total MI for each gene pair
teResSum = sum(teRes, dims=3)[:,:,1]
#zero the diagonal
for i=1:10
    for j=1:10
        if i==j
            teResSum[i,j] = 0
        end
    end
end

rankVecTE = Vector{Vector{Int}}(undef, 10)
for i=1:10
    causalGenes = findall(x->x != 0, geneTFAMat[:,i])
    ranks = findall(in(causalGenes), reverse(sortperm(teResSum[:,i])))
    rankVecTE[i] = ranks
end


using JLD2
mcdsLarge = load_object("Results/multipleCellSimulation_100Genes.jld2")
