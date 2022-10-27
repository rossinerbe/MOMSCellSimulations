include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

#change splicing rate
splicingGS = Vector{Gene}(undef, length(testGS))
for i=1:5
    splicingGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        1, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=6:10
    splicingGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        10, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end

sIC = generateInitialConditions(splicingGS)
testSimSplicing = cellSimulation(splicingGS, sIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimSplicing.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimSplicing.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimSplicing.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimSplicing.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp1_ProteinExpression.pdf")



splicingGS = Vector{Gene}(undef, length(testGS))
for i=1:5
    splicingGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        1, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end
for i=6:10
    splicingGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        25, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end

sIC = generateInitialConditions(splicingGS)
testSimSplicing = cellSimulation(splicingGS, sIC, 1000)


plot(x, transpose(testSimSplicing.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp2_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimSplicing.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp2_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimSplicing.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp2_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimSplicing.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testSp2_ProteinExpression.pdf")
