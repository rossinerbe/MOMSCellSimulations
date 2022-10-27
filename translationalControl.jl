include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

#remove translational inhibitors

tiGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    tiGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, [],
        testGS[i].ProteinDegradationFactors)
end

tiIC = generateInitialConditions(tiGS)
testSimTI = cellSimulation(tiGS, tiIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimTI.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTI1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimTI.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTI1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTI.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTI1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTI.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTI1_ProteinExpression.pdf")


#equalize TIP
tipGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    tipGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        0.4, testGS[i].TranslationalInhibitor,
        testGS[i].ProteinDegradationFactors)
end

tipIC = generateInitialConditions(tipGS)
testSimTIP = cellSimulation(tipGS, tipIC, 1000)


plot(x, transpose(testSimTIP.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIP1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimTIP.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIP1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTIP.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIP1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTIP.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIP1_ProteinExpression.pdf")


#both
tiptiGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    tiptiGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        0.4, [], testGS[i].ProteinDegradationFactors)
end

tiptiIC = generateInitialConditions(tiptiGS)
testSimTIPTI = cellSimulation(tiptiGS, tiptiIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimTIPTI.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIPTI1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimTIPTI.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIPTI1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTIPTI.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIPTI1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimTIPTI.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testTIPTI1_ProteinExpression.pdf")


#remove protein degradation factors
pdfGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    pdfGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
        testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
        testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
        testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
        testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
        [])
end

pdfIC = generateInitialConditions(pdfGS)
testSimPDF = cellSimulation(pdfGS, pdfIC, 1000)


plot(x, transpose(testSimPDF.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPDF1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimPDF.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPDF1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimPDF.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPDF1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimPDF.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testPDF1_ProteinExpression.pdf")
