include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

btpGS = Vector{Gene}(undef, length(testGS))
for i=1:5
    btpGS[i] = Gene(testGS[i].Name, 0.1,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=6:10
    btpGS[i] = Gene(testGS[i].Name, 0.005,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end


btpIC = generateInitialConditions(btpGS)
testSimBTP = cellSimulation(btpGS, btpIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimBTP.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimBTP.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimBTP.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimBTP.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP1_ProteinExpression.pdf")



btpGS = Vector{Gene}(undef, length(testGS))
for i=1
    btpGS[i] = Gene(testGS[i].Name, 0.1,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=2:10
    btpGS[i] = Gene(testGS[i].Name, 0.005,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end


btpIC = generateInitialConditions(btpGS)
testSimBTP = cellSimulation(btpGS, btpIC, 1000)


plot(x, transpose(testSimBTP.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP2_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimBTP.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP2_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimBTP.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP2_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimBTP.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testBTP2_ProteinExpression.pdf")
