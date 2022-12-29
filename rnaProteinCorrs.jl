include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true, seed=1)
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

#correlation of RNA and protein levels
using Statistics
RPCorrs = Vector{Float64}(undef, length(testGS))
for i=1:length(RPCorrs)
    RPCorrs[i] = cor(testSim.SplicedRNA[i,:], testSim.ProteinExpression[i,:])
end

#compute cross (time lagged) correlation
using StatsBase
RPCrossCorrs = Vector{Float64}(undef, length(testGS))
for i=1:length(RPCrossCorrs)
    RPCrossCorrs[i] = crosscor(float.(testSim.SplicedRNA[i,:]), float.(testSim.ProteinExpression[i,:]),[(testGS[i].SplicingRate+3)])[1]
end

#get correlations of RNA and protein for 100 different random gene sets
corrRes = Array{Float64}(undef, 10, 100)
for i=1:100
    testGS = generateGeneSet(10, allowFeedbackLoops = true, seed=i)
    testIC = generateInitialConditions(testGS)
    testSim = cellSimulation(testGS, testIC, 1000)
    for j=1:10
        corrRes[j,i] = cor(testSim.SplicedRNA[j,:], testSim.ProteinExpression[j,:])
    end
end

mean(corrRes)
median(corrRes)
length(findall(x->x < 0, corrRes))
nquantile(vec(corrRes), 5)
minimum(corrRes)
maximum(corrRes)

histogram(sort(vec(corrRes)), legend = false)
savefig("C:/Users/rossi/Pictures/CellSimulations/RNAProteinCorrsHist.pdf")

simMeds = mean(corrRes, dims =1)
histogram(sort(vec(simMeds)), legend = false)
savefig("C:/Users/rossi/Pictures/CellSimulations/RNAProteinCorrsSimMediansHist.pdf")

simMeans = mean(corrRes, dims =1)
histogram(sort(vec(simMeans)), legend = false)
savefig("C:/Users/rossi/Pictures/CellSimulations/RNAProteinCorrsSimMeansHist.pdf")

#least correlated simulation plots
findall(x->x == minimum(simMeans), simMeans)
testGS = generateGeneSet(10, allowFeedbackLoops = true, seed=70)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

x= collect(1:1000)
plot(x, transpose(testSim.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_lowestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_lowestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_lowestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_lowestRNAProteinCorr.pdf")


#most correlated simulation plots
findall(x->x == maximum(simMeans), simMeans)
testGS = generateGeneSet(10, allowFeedbackLoops = true, seed=96)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

x= collect(1:1000)
plot(x, transpose(testSim.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_highestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_highestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_highestRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_highestRNAProteinCorr.pdf")

#median level of correlation
findmin(abs.(simMeans.-median(simMeans)))[2]
testGS = generateGeneSet(10, allowFeedbackLoops = true, seed=46)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

plot(x, transpose(testSim.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_medianRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_medianRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_medianRNAProteinCorr.pdf")

plot(x, log1p.(transpose(testSim.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_medianRNAProteinCorr.pdf")


#compare to information theory - transfer entropy
using TransferEntropy
transferentropy(testSim.SplicedRNA[1,:], testSim.ProteinExpression[1,:], VisitationFrequency(RectangularBinning(250)))
