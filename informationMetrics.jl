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

#create data set with many cells from different gene sets
#made using different random seeds
mgsData = Vector{Any}(undef, 100)
for i=1:100
    tmpGS = generateGeneSet(10, seed = i)
    tmpIC = generateInitialConditions(tmpGS)
    mgsData[i] = cellSimulation(tmpGS, tmpIC, 500)
end

#and with a larger number of genes
mgsData100Genes = Vector{Any}(undef, 100)
for i=1:100
    tmpGS = generateGeneSet(100, seed = i)
    tmpIC = generateInitialConditions(tmpGS)
    mgsData100Genes[i] = cellSimulation(tmpGS, tmpIC, 500)
end

#find mutual information rna-rna, rna-protein, protein-protein
est =VisitationFrequency(RectangularBinning(200))
#for simulated single cell data set
miSCDS = (RNA_RNA=zeros(10, 10), RNA_Protein=zeros(10, 10),
 Protein_Protein=zeros(10, 10))
for i=1:10
    for j=1:10
            miSCDS.RNA_RNA[i,j] = mutualinfo(scDS.SplicedRNA[i,:], scDS.SplicedRNA[j,:], est)
            miSCDS.RNA_Protein[i,j] = mutualinfo(scDS.SplicedRNA[i,:], scDS.ProteinExpression[j,:], est)
            miSCDS.Protein_Protein[i,j] = mutualinfo(scDS.ProteinExpression[i,:], scDS.ProteinExpression[j,:], est)
    end
end

#save result
using JLD2
save_object("Results/mutualInfo_MultiCellSnapshots.jld2", miSCDS)

#compare to RNA protein correlation
rpCorr = zeros(10,10)
for i=1:10
    for j=1:10
            rpCorr[i,j] = cor(scDS.SplicedRNA[i,:], scDS.ProteinExpression[j,:])
    end
end

using LinearAlgebra
median(diag(rpCorr))
mean(diag(rpCorr))

#save result
using JLD2
save_object("Results/rnaProteinCorr_MultiCellSnapshots.jld2", rpCorr)

#larger data set
miSCDSLarge = (RNA_RNA=zeros(100, 100), RNA_Protein=zeros(100, 100),
 Protein_Protein=zeros(100, 100))
for i=1:100
    for j=1:100
            miSCDSLarge.RNA_RNA[i,j] = mutualinfo(scDSLarge.SplicedRNA[i,:], scDSLarge.SplicedRNA[j,:], est)
            miSCDSLarge.RNA_Protein[i,j] = mutualinfo(scDSLarge.SplicedRNA[i,:], scDSLarge.ProteinExpression[j,:], est)
            miSCDSLarge.Protein_Protein[i,j] = mutualinfo(scDSLarge.ProteinExpression[i,:], scDSLarge.ProteinExpression[j,:], est)
    end
end

#save result

save_object("Results/mutualInfo_MultiCellSnapshots_100Genes.jld2", miSCDSLarge)

#compare to RNA protein correlation
rpCorrLarge = zeros(100,100)
for i=1:100
    for j=1:100
            rpCorrLarge[i,j] = cor(scDSLarge.SplicedRNA[i,:], scDSLarge.ProteinExpression[j,:])
    end
end


rnaCorrProtCorrs = diag(rpCorrLarge)
median(rnaCorrProtCorrs)
mean(rnaCorrProtCorrs)

#save result
save_object("Results/rnaProteinCorr_MultiCellSnapshots_100Genes.jld2", rpCorrLarge)

#look at time matched data
miMGS = Vector{Any}(undef, 100)
for k=1:100
    tmp = (RNA_RNA=zeros(10, 10), RNA_Protein=zeros(10, 10),
     Protein_Protein=zeros(10, 10))
    for i=1:10
        for j=1:10
                tmp.RNA_RNA[i,j] = mutualinfo(mgsData[k].SplicedRNA[i,:], mgsData[k].SplicedRNA[j,:], est)
                tmp.RNA_Protein[i,j] = mutualinfo(mgsData[k].SplicedRNA[i,:], mgsData[k].ProteinExpression[j,:], est)
                tmp.Protein_Protein[i,j] = mutualinfo(mgsData[k].ProteinExpression[i,:], mgsData[k].ProteinExpression[j,:], est)
        end
    end
    miMGS[k] = tmp
end

miMGS100Genes = Vector{Any}(undef, 100)
for k=1:100
    tmp = (RNA_RNA=zeros(100, 100), RNA_Protein=zeros(100, 100),
     Protein_Protein=zeros(100, 100))
    for i=1:100
        for j=1:100
            try
                tmp.RNA_RNA[i,j] = mutualinfo(mgsData100Genes[k].SplicedRNA[i,:], mgsData100Genes[k].SplicedRNA[j,:], est)
                tmp.RNA_Protein[i,j] = mutualinfo(mgsData100Genes[k].SplicedRNA[i,:], mgsData100Genes[k].ProteinExpression[j,:], est)
                tmp.Protein_Protein[i,j] = mutualinfo(mgsData100Genes[k].ProteinExpression[i,:], mgsData100Genes[k].ProteinExpression[j,:], est)
            catch
                tmp.RNA_RNA[i,j] = 0
                tmp.RNA_Protein[i,j] = 0
                tmp.Protein_Protein[i,j] = 0
            end
        end
    end
    miMGS100Genes[k] = tmp
end

#for each rna-protein relationship, calculate the proportion of genes where there
#is another gene with more MI with a protein than the RNA that is used to produce it
higherMICount = 0
for i=1:10
    if length(findall(x->x > miSCDS.RNA_Protein[i,i], miSCDS.RNA_Protein[:,i])) > 0
        higherMICount+=1
    end
end

higherGeneRecord = Vector{Any}(undef,100)
higherMICountL = 0
for i=1:100
    if length(findall(x->x > miSCDSLarge.RNA_Protein[i,i], miSCDSLarge.RNA_Protein[:,i])) > 0
        higherMICountL+=1
        higherGeneRecord[i] = findall(x->x > miSCDSLarge.RNA_Protein[i,i], miSCDSLarge.RNA_Protein[:,i])
    else
        higherGeneRecord[i] = []
    end
end


#find transfer entropy rna->rna, rna->protein, protein->protein, protein->rna
teSCDS = (RNA_RNA=zeros(10, 10), RNA_Protein=zeros(10, 10),
 Protein_Protein=zeros(10, 10), Protein_RNA = zeros(10, 10))
for i=1:10
    for j=1:10
            teSCDS.RNA_RNA[i,j] = transferentropy(scDS.SplicedRNA[i,:], scDS.SplicedRNA[j,:], est)
            teSCDS.RNA_Protein[i,j] = transferentropy(scDS.SplicedRNA[i,:], scDS.ProteinExpression[j,:], est)
            teSCDS.Protein_Protein[i,j] = transferentropy(scDS.ProteinExpression[i,:], scDS.ProteinExpression[j,:], est)
            teSCDS.Protein_RNA[i,j] = transferentropy(scDS.ProteinExpression[i,:], scDS.SplicedRNA[j,:], est)
    end
end

save_object("Results/transferEntropy_MultiCellSnapshots.jld2", teSCDS)

#larger data set
teSCDSLarge = (RNA_RNA=zeros(100, 100), RNA_Protein=zeros(100, 100),
 Protein_Protein=zeros(100, 100), Protein_RNA = zeros(100, 100))
for i=1:100
    for j=1:100
            teSCDSLarge.RNA_RNA[i,j] = transferentropy(scDSLarge.SplicedRNA[i,:], scDSLarge.SplicedRNA[j,:], est)
            teSCDSLarge.RNA_Protein[i,j] = transferentropy(scDSLarge.SplicedRNA[i,:], scDSLarge.ProteinExpression[j,:], est)
            teSCDSLarge.Protein_Protein[i,j] = transferentropy(scDSLarge.ProteinExpression[i,:], scDSLarge.ProteinExpression[j,:], est)
            teSCDSLarge.Protein_RNA[i,j] = transferentropy(scDSLarge.ProteinExpression[i,:], scDSLarge.SplicedRNA[j,:], est)
    end
end

#save result
save_object("Results/transferEntropy_MultiCellSnapshots_100Genes.jld2", teSCDSLarge)

#compare to adjacency matrix, mainly of transcriptional regulators for gene regulatory
#network inference of TF-gene relationships
geneTFAMat = getTFGeneAMat(testGS)
rankVec = Vector{Vector{Int}}(undef, 10)
for i=1:10
    causalGenes = findall(x->x != 0, geneTFAMat[:,i])
    #set diagonal of mi matrix to zero
    diag(miSCDS.RNA_RNA) .= 0
    ranks = findall(in(causalGenes), reverse(sortperm(miSCDS.RNA_RNA[:,i])))
    rankVec[i] = ranks
end


geneTFAMatLarge = getTFGeneAMat(testGSLarge)
rankVecLarge = Vector{Vector{Int}}(undef, 100)
for i=1:100
    causalGenes = findall(x->x != 0, geneTFAMatLarge[:,i])
    #set diagonal of mi matrix to zero
    diag(miSCDSLarge.RNA_RNA) .= 0
    ranks = findall(in(causalGenes), reverse(sortperm(miSCDSLarge.RNA_RNA[:,i])))
    rankVecLarge[i] = ranks
end

#bar plots of how many genes are in each rank
ranks = vcat(rankVec...)
rankCounts = countmap(ranks)
rankCountsOrdered = Vector{Int}(undef, 10)
for i=1:10
    ind = findall(x->x == i, collect(keys(rankCounts)))
    if length(ind) == 0
        rankCountsOrdered[i] = 0
    else
        rankCountsOrdered[i] = collect(values(rankCounts))[ind[1]]
    end
end

ranksLarge = vcat(rankVecLarge...)
rankCountsLarge = countmap(ranksLarge)
rankCountsOrderedL = Vector{Int}(undef, 100)
for i=1:100
    ind = findall(x->x == i, collect(keys(rankCountsLarge)))
    if length(ind) == 0
        rankCountsOrderedL[i] = 0
    else
        rankCountsOrderedL[i] = collect(values(rankCountsLarge))[ind[1]]
    end
end

using Plots
bar(rankCountsOrdered, orientation=:h, yticks=(1:10), yflip=true,
 xlabel = "Number of Genes", ylabel = "Mutual information rank", legend = false)
savefig("MutualInformation/MIRankofDirectTranscriptionalRegulators_10gene.pdf")
bar(rankCountsOrderedL, orientation=:h, yticks=(1:100), yflip=true,
 xlabel = "Number of Genes", ylabel = "Mutual information rank", legend = false)
savefig("MutualInformation/MIRankofDirectTranscriptionalRegulators_100gene.pdf")


#bar plots for TE
rankVecTE = Vector{Vector{Int}}(undef, 10)
for i=1:10
    causalGenes = findall(x->x != 0, geneTFAMat[:,i])
    #set diagonal of te matrix to zero
    diag(teSCDS.RNA_RNA) .= 0
    ranks = findall(in(causalGenes), reverse(sortperm(teSCDS.RNA_RNA[:,i])))
    rankVecTE[i] = ranks
end


rankVecTELarge = Vector{Vector{Int}}(undef, 100)
for i=1:100
    causalGenes = findall(x->x != 0, geneTFAMatLarge[:,i])
    #set diagonal of mi matrix to zero
    diag(teSCDSLarge.RNA_RNA) .= 0
    ranks = findall(in(causalGenes), reverse(sortperm(teSCDSLarge.RNA_RNA[:,i])))
    rankVecTELarge[i] = ranks
end

#bar plots of how many genes are in each rank
ranksTE = vcat(rankVecTE...)
rankCountsTE = countmap(ranksTE)
rankCountsOrderedTE = Vector{Int}(undef, 10)
for i=1:10
    ind = findall(x->x == i, collect(keys(rankCountsTE)))
    if length(ind) == 0
        rankCountsOrderedTE[i] = 0
    else
        rankCountsOrderedTE[i] = collect(values(rankCountsTE))[ind[1]]
    end
end

ranksTELarge = vcat(rankVecTELarge...)
rankCountsTELarge = countmap(ranksTELarge)
rankCountsOrderedTEL = Vector{Int}(undef, 100)
for i=1:100
    ind = findall(x->x == i, collect(keys(rankCountsTELarge)))
    if length(ind) == 0
        rankCountsOrderedTEL[i] = 0
    else
        rankCountsOrderedTEL[i] = collect(values(rankCountsTELarge))[ind[1]]
    end
end

bar(rankCountsOrderedTE, orientation=:h, yticks=(1:10), yflip=true,
 xlabel = "Number of Genes", ylabel = "Transfer entropy rank", legend = false)
savefig("TransferEntropy/TERankofDirectTranscriptionalRegulators_10gene.pdf")
bar(rankCountsOrderedTEL, orientation=:h, yticks=(1:100), yflip=true,
 xlabel = "Number of Genes", ylabel = "Transfer entropy rank", legend = false)
savefig("TransferEntropy/TERankofDirectTranscriptionalRegulators_100gene.pdf")
