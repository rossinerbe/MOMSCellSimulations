include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

#change protein degradation rate

phlGS = Vector{Gene}(undef, length(testGS))
for i=1:5
    phlGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, 50,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=6:10
    phlGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, 1500,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end


phlIC = generateInitialConditions(phlGS)
testSimPHL = cellSimulation(phlGS, phlIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimPHL.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPHL1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimPHL.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPHL1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimPHL.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPHLP1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimPHL.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPHL1_ProteinExpression.pdf")


#change RNA degradation rate
rhlGS = Vector{Gene}(undef, length(testGS))
for i=1:5
    rhlGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, 1200, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=6:10
    rhlGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, 40, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end


rhlIC = generateInitialConditions(rhlGS)
testSimRHL = cellSimulation(rhlGS, rhlIC, 1000)


plot(x, transpose(testSimRHL.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testRHL1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimRHL.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testRHL1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimRHL.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testRHL1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimRHL.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testRHL1_ProteinExpression.pdf")
