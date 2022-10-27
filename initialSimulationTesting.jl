include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

using Plots
x= collect(1:1000)
plot(x, transpose(testSim.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState.pdf")

plot(x, log1p.(transpose(testSim.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts.pdf")

plot(x, log1p.(transpose(testSim.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSim.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression.pdf")

#what does it look like with a different random seed?
testSim2 = cellSimulation(testGS, testIC, 1000, seed = 42)

x= collect(1:1000)
plot(x, transpose(testSim2.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_seed2.pdf")

plot(x, log1p.(transpose(testSim2.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_seed2.pdf")

plot(x, log1p.(transpose(testSim2.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_seed2.pdf")

plot(x, log1p.(transpose(testSim2.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_seed2.pdf")


#how much effect does simulation seed have?
#how much effect does changing ICS have?
#how different is another random gene set with only seed changed in generation?
#how about if you change different parameters?
#how to determine size of effect?

#create the data
simSeedData = Vector{Any}(undef, 100)
for i=1:100
    simSeedData[i] = cellSimulation(testGS, testIC, 1000, seed = i)
end

simICSData = Vector{Any}(undef, 100)
for i=1:100
    tmpIC = generateInitialConditions(testGS, seed = i)
    simICSData[i] = cellSimulation(testGS, tmpIC, 1000)
end

simGSData = Vector{Any}(undef, 100)
for i=1:100
    tmpGS = generateGeneSet(10, seed = i)
    tmpIC = generateInitialConditions(tmpGS)
    simGSData[i] = cellSimulation(tmpGS, tmpIC, 1000)
end

using JLD2
save_object("D:/FertigLabLargeData/CellSimulations/Results/simulatedSeedData.jld2", simSeedData)
save_object("D:/FertigLabLargeData/CellSimulations/Results/simulatedICSData.jld2", simICSData)
save_object("D:/FertigLabLargeData/CellSimulations/Results/simulatedGSData.jld2", simGSData)


#find the correlation for each modality
using Statistics, LinearAlgebra
simSeedCorrs = Vector{Any}(undef, length(simSeedData))
for i=1:length(simSeedData)
    tmpSplicedCor = diag(cor(simSeedData[i].SplicedRNA, testSim.SplicedRNA, dims = 2))
    tmpUnsplicedCor = diag(cor(simSeedData[i].UnsplicedRNA, testSim.UnsplicedRNA, dims = 2))
    tmpProteinCor = diag(cor(simSeedData[i].ProteinExpression, testSim.ProteinExpression, dims = 2))
    tmpEpigeneticCor = diag(cor(simSeedData[i].EpigeneticState, testSim.EpigeneticState, dims = 2))
    simSeedCorrs[i] = (SplicedRNA = tmpSplicedCor, UnsplicedRNA = tmpUnsplicedCor,
    Protein = tmpProteinCor, Epigenetic = tmpEpigeneticCor)
end

#changing the random seed can lead to different trajectories for some genes. Most genes follow
#similar trajectories, but usually 1-2 diverge noticeably

icsCorrs = Vector{Any}(undef, length(simICSData))
for i=1:length(simICSData)
    tmpSplicedCor = diag(cor(simICSData[i].SplicedRNA, testSim.SplicedRNA, dims = 2))
    tmpUnsplicedCor = diag(cor(simICSData[i].UnsplicedRNA, testSim.UnsplicedRNA, dims = 2))
    tmpProteinCor = diag(cor(simICSData[i].ProteinExpression, testSim.ProteinExpression, dims = 2))
    tmpEpigeneticCor = diag(cor(simICSData[i].EpigeneticState, testSim.EpigeneticState, dims = 2))
    icsCorrs[i] = (SplicedRNA = tmpSplicedCor, UnsplicedRNA = tmpUnsplicedCor,
    Protein = tmpProteinCor, Epigenetic = tmpEpigeneticCor)
end

#plot an example
plot(x, transpose(simICSData[1].EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_diffICS.pdf")

plot(x, log1p.(transpose(simICSData[1].SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_diffICS.pdf")

plot(x, log1p.(transpose(simICSData[1].UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_diffICS.pdf")

plot(x, log1p.(transpose(simICSData[1].ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_diffICS.pdf")

#second example
plot(x, transpose(simICSData[22].EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_diffICS2.pdf")

plot(x, log1p.(transpose(simICSData[22].SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_diffICS2.pdf")

plot(x, log1p.(transpose(simICSData[22].UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_diffICS2.pdf")

plot(x, log1p.(transpose(simICSData[22].ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_diffICS2.pdf")

#correlations is much less for different ICS, but the shape of the overall progression through
#time is still consistent for most genes, though it depends on the example


#gene set
gsCorrs = Vector{Any}(undef, length(simGSData))
for i=1:length(simGSData)
    tmpSplicedCor = diag(cor(simGSData[i].SplicedRNA, testSim.SplicedRNA, dims = 2))
    tmpUnsplicedCor = diag(cor(simGSData[i].UnsplicedRNA, testSim.UnsplicedRNA, dims = 2))
    tmpProteinCor = diag(cor(simGSData[i].ProteinExpression, testSim.ProteinExpression, dims = 2))
    tmpEpigeneticCor = diag(cor(simGSData[i].EpigeneticState, testSim.EpigeneticState, dims = 2))
    gsCorrs[i] = (SplicedRNA = tmpSplicedCor, UnsplicedRNA = tmpUnsplicedCor,
    Protein = tmpProteinCor, Epigenetic = tmpEpigeneticCor)
end


#plot example
plot(x, transpose(simGSData[1].EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_diffGS.pdf")

plot(x, log1p.(transpose(simGSData[1].SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_diffGS.pdf")

plot(x, log1p.(transpose(simGSData[1].UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_diffGS.pdf")

plot(x, log1p.(transpose(simGSData[1].ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_diffGS.pdf")

#second example
plot(x, transpose(simGSData[99].EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_diffGS2.pdf")

plot(x, log1p.(transpose(simGSData[99].SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_diffGS2.pdf")

plot(x, log1p.(transpose(simGSData[99].UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_diffGS2.pdf")

plot(x, log1p.(transpose(simGSData[99].ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_diffGS2.pdf")


##100 genes
testGS = generateGeneSet(100, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

plot(x, transpose(testSim.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testEpigeneticState_100Genes.pdf")

plot(x, log1p.(transpose(testSim.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSplicedCounts_100Genes.pdf")

plot(x, log1p.(transpose(testSim.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testUnsplicedCounts_100Genes.pdf")

plot(x, log1p.(transpose(testSim.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testProteinExpression_100Genes.pdf")
