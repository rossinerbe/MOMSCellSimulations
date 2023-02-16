include("simulateCells.jl");

using TransferEntropy
using Statistics

#make 100  10 gene data sets over 500 time points with one hundred different random
#regulatory interaction networks
mgsData = Vector{Any}(undef, 100)
for i=1:100
    tmpGS = generateGeneSet(10, seed = i)
    tmpIC = generateInitialConditions(tmpGS)
    mgsData[i] = cellSimulation(tmpGS, tmpIC, 500)
end

#and with a 100 genes
mgsData100Genes = Vector{Any}(undef, 100)
for i=1:100
    tmpGS = generateGeneSet(100, seed = i)
    tmpIC = generateInitialConditions(tmpGS)
    mgsData100Genes[i] = cellSimulation(tmpGS, tmpIC, 500)
end

#MI estimator
est =VisitationFrequency(RectangularBinning(200))

#MI beetween epigenetic state and RNA, epg state and protein, and rna and protein
#for ten gene data
miMGS = Vector{Any}(undef, 100)
for k=1:100
    tmp = (EPG_RNA=zeros(10, 10), EPG_Protein=zeros(10, 10),
     RNA_Protein=zeros(10, 10))
    for i=1:10
        for j=1:10
            if var(mgsData[k].EpigeneticState[i,:]) == 0.0
                tmp.EPG_RNA[i,j] = 0.0
                tmp.EPG_Protein[i,j] = 0.0
                tmp.RNA_Protein[i,j] = mutualinfo(mgsData[k].SplicedRNA[i,:], mgsData[k].ProteinExpression[j,:], est)
            else
                tmp.EPG_RNA[i,j] = mutualinfo(mgsData[k].EpigeneticState[i,:], mgsData[k].SplicedRNA[j,:], est)
                tmp.EPG_Protein[i,j] = mutualinfo(mgsData[k].EpigeneticState[i,:], mgsData[k].ProteinExpression[j,:], est)
                tmp.RNA_Protein[i,j] = mutualinfo(mgsData[k].SplicedRNA[i,:], mgsData[k].ProteinExpression[j,:], est)
            end
        end
    end
    miMGS[k] = tmp
end


#mean and median MI from EPG to RNA and protein vs RNA-protein
#from each sim
means = Matrix{Float64}(undef, 100,3)
medians = Matrix{Float64}(undef, 100,3)

#don't use zeroes because the point is to look at impact of chromatin variation
#and those are the zero variation genes
for i=1:100
    means[i,1] = mean(miMGS[i].EPG_RNA[findall(x->x > 0, miMGS[i].EPG_RNA)])
    means[i,2] = mean(miMGS[i].EPG_Protein[findall(x->x > 0, miMGS[i].EPG_Protein)])
    means[i,3] = mean(miMGS[i].RNA_Protein[findall(x->x > 0, miMGS[i].RNA_Protein)])

    medians[i,1] = median(miMGS[i].EPG_RNA[findall(x->x > 0, miMGS[i].EPG_RNA)])
    medians[i,2] = median(miMGS[i].EPG_Protein[findall(x->x > 0, miMGS[i].EPG_Protein)])
    medians[i,3] = median(miMGS[i].RNA_Protein[findall(x->x > 0, miMGS[i].RNA_Protein)])
end

#visualize
using Plots
using StatsPlots
boxplot(["EPG_RNA" "EPG_Protein" "RNA_Protein"], means, legend = false,
    ylabel = "Mean Mutual Information")
savefig("MutualInformation/MeanMutualInformationEpgvAllBoxplot.pdf")

boxplot(["EPG_RNA" "EPG_Protein" "RNA_Protein"], medians, legend = false,
    ylabel = "Median Mutual Information")
savefig("MutualInformation/MedianMutualInformationEpgvAllBoxplot.pdf")


#repeat for 100 gene sims
miMGS100Genes = Vector{Any}(undef, 100)
for k=1:100
    tmp = (EPG_RNA=zeros(100, 100), EPG_Protein=zeros(100, 100),
     RNA_Protein=zeros(100, 100))
    for i=1:100
        for j=1:100
            if var(mgsData100Genes[k].EpigeneticState[i,:]) == 0.0
                tmp.EPG_RNA[i,j] = 0.0
                tmp.EPG_Protein[i,j] = 0.0
                if var(mgsData100Genes[k].SplicedRNA[j,:]) == 0.0 || var(mgsData100Genes[k].SplicedRNA[i,:]) == 0.0
                    tmp.RNA_Protein[i,j] = 0.0
                else
                    tmp.RNA_Protein[i,j] = mutualinfo(mgsData100Genes[k].SplicedRNA[i,:], mgsData100Genes[k].ProteinExpression[j,:], est)
                end
            elseif var(mgsData100Genes[k].SplicedRNA[j,:]) == 0.0 || var(mgsData100Genes[k].SplicedRNA[i,:]) == 0.0
                tmp.EPG_RNA[i,j] = 0.0
                tmp.EPG_Protein[i,j] = 0.0
                tmp.RNA_Protein[i,j] = 0.0
            elseif var(mgsData100Genes[k].ProteinExpression[j,:]) == 0.0
                tmp.EPG_RNA[i,j] = 0.0
                tmp.EPG_Protein[i,j] = 0.0
                tmp.RNA_Protein[i,j] = 0.0
            else
                tmp.EPG_RNA[i,j] = mutualinfo(mgsData100Genes[k].EpigeneticState[i,:], mgsData100Genes[k].SplicedRNA[j,:], est)
                tmp.EPG_Protein[i,j] = mutualinfo(mgsData100Genes[k].EpigeneticState[i,:], mgsData100Genes[k].ProteinExpression[j,:], est)
                tmp.RNA_Protein[i,j] = mutualinfo(mgsData100Genes[k].SplicedRNA[i,:], mgsData100Genes[k].ProteinExpression[j,:], est)
            end
        end
    end
    miMGS100Genes[k] = tmp
end


#mean and median MI from EPG to RNA and protein vs RNA-protein
#from each sim
means100Genes = Matrix{Float64}(undef, 100,3)
medians100Genes = Matrix{Float64}(undef, 100,3)

#don't use zeroes because the point is to look at impact of chromatin variation
#and those are the zero variation genes
for i=1:100
    means100Genes[i,1] = mean(miMGS100Genes[i].EPG_RNA[findall(x->x > 0, miMGS100Genes[i].EPG_RNA)])
    means100Genes[i,2] = mean(miMGS100Genes[i].EPG_Protein[findall(x->x > 0, miMGS100Genes[i].EPG_Protein)])
    means100Genes[i,3] = mean(miMGS100Genes[i].RNA_Protein[findall(x->x > 0, miMGS100Genes[i].RNA_Protein)])

    medians100Genes[i,1] = median(miMGS100Genes[i].EPG_RNA[findall(x->x > 0, miMGS100Genes[i].EPG_RNA)])
    medians100Genes[i,2] = median(miMGS100Genes[i].EPG_Protein[findall(x->x > 0, miMGS100Genes[i].EPG_Protein)])
    medians100Genes[i,3] = median(miMGS100Genes[i].RNA_Protein[findall(x->x > 0, miMGS100Genes[i].RNA_Protein)])
end

#visualize
boxplot(["EPG_RNA" "EPG_Protein" "RNA_Protein"], means100Genes, legend = false,
    ylabel = "Mean Mutual Information")
savefig("MutualInformation/MeanMutualInformationEpgvAllBoxplot_100Genes.pdf")

boxplot(["EPG_RNA" "EPG_Protein" "RNA_Protein"], medians100Genes, legend = false,
    ylabel = "Median Mutual Information")
savefig("MutualInformation/MedianMutualInformationEpgvAllBoxplot_100Genes.pdf")

#how well does epigenetic state predict RNA level in a simulated single cell data set?
#specifically, if accessibility is 0, is gene bottom third
#if 1, mid third, and 2 top 3rd
#create data set with many cells from the same gene set
testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
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

epgGeneVar = var(scDS.EpigeneticState, dims=2)
sRNAGeneVar = var(scDS.SplicedRNA, dims=2)
proteinGeneVar = var(scDS.ProteinExpression, dims=2)
cor(epgGeneVar, sRNAGeneVar)
cor(epgGeneVar, proteinGeneVar)
