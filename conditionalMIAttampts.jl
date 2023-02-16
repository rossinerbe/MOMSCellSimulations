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

#select random snapshots of cells to use as "realistic" single cell data set
#don't use first 10 time points when simulation is just starting because it is
#less stable there and a worse approximation of normal cell dynamics
#Random and StatsBase should already be loaded from simulateCells.jl
Random.seed!(123)
scDS = (EpigeneticState = Matrix{Int}(undef, 10, 2000),
        SplicedRNA = Matrix{Int}(undef, 10, 2000),
        ProteinExpression = Matrix{Int}(undef, 10, 2000))
for i=1:2000
    tp = sample(10:500)
    scDS.EpigeneticState[:,i] = mcds[i].EpigeneticState[:,tp]
    scDS.SplicedRNA[:,i] = mcds[i].SplicedRNA[:,tp]
    scDS.ProteinExpression[:,i] = mcds[i].ProteinExpression[:,tp]
end

#one with 100 genes
testGSLarge = generateGeneSet(100, allowFeedbackLoops = true)
testICLarge = generateInitialConditions(testGSLarge)
mcdsLarge = makeMultiCellDataSet(testGSLarge, testICLarge, 500, 2000)
save_object("Results/multipleCellSimulation_100Genes.jld2", mcdsLarge)
Random.seed!(123)
scDSLarge = (EpigeneticState = Matrix{Int}(undef, 100, 2000),
        SplicedRNA = Matrix{Int}(undef, 100, 2000),
        ProteinExpression = Matrix{Int}(undef, 100, 2000))
for i=1:2000
    tp = sample(10:500)
    scDSLarge.EpigeneticState[:,i] = mcdsLarge[i].EpigeneticState[:,tp]
    scDSLarge.SplicedRNA[:,i] = mcdsLarge[i].SplicedRNA[:,tp]
    scDSLarge.ProteinExpression[:,i] = mcdsLarge[i].ProteinExpression[:,tp]
end


est =VisitationFrequency(RectangularBinning(200))
miSCDS = (RNA_RNA_ATAC=zeros(10, 10), RNA_RNA_Protein=zeros(10, 10))
for i=1:10
    for j=1:10
            miSCDS.RNA_RNA_ATAC[i,j] = conditional_mutualinfo(scDS.SplicedRNA[i,:], scDS.SplicedRNA[j,:], scDS.EpigeneticState[i,:], est)
            miSCDS.RNA_Protein[i,j] = mutualinfo(scDS.SplicedRNA[i,:], scDS.ProteinExpression[j,:], est)
            miSCDS.Protein_Protein[i,j] = mutualinfo(scDS.ProteinExpression[i,:], scDS.ProteinExpression[j,:], est)
    end
end
