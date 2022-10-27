include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

#change transcriptional network to introduce negative feedback on genes with highest protein expression

nfGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    if i == 1 || i==2 || i==7
        gn = "Gene"*string(i)
        nfGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
            testGS[i].TranscriptionActivators, vcat(testGS[i].TranscriptionRepressors, gn),
            testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
            testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
            testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
            testGS[i].ProteinDegradationFactors)
    else
        nfGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
            testGS[i].TranscriptionActivators, testGS[i].TranscriptionRepressors,
            testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
            testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
            testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
            testGS[i].ProteinDegradationFactors)
        end
end



nfIC = generateInitialConditions(nfGS)
testSimNF = cellSimulation(nfGS,nfIC, 1000)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimNF.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNF1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimNF.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNF1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimNF.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNF1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimNF.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNF1_ProteinExpression.pdf")


#try positive feedback for lowly expressed genes, in tandem with negative feedback for up-expressed genes

npfGS = Vector{Gene}(undef, length(testGS))
for i=1:10
    gn = "Gene"*string(i)
    if i == 1 || i==2 || i==7
        npfGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
            testGS[i].TranscriptionActivators, vcat(testGS[i].TranscriptionRepressors, gn),
            testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
            testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
            testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
            testGS[i].ProteinDegradationFactors)
    else
        npfGS[i] = Gene(testGS[i].Name, testGS[i].BaseTranscriptionProbability,
            vcat(testGS[i].TranscriptionActivators, gn), testGS[i].TranscriptionRepressors,
            testGS[i].EpigeneticUpregulators, testGS[i].EpigeneticDownregulators,
            testGS[i].SplicingRate, testGS[i].RNAHalfLife, testGS[i].ProteinHalfLife,
            testGS[i].TranslationInitiationProbability, testGS[i].TranslationalInhibitor,
            testGS[i].ProteinDegradationFactors)
        end
end


npfIC = generateInitialConditions(npfGS)
testSimNPF = cellSimulation(npfGS,npfIC, 1000)

plot(x, transpose(testSimNPF.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimNPF.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimNPF.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimNPF.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF1_ProteinExpression.pdf")


#longer run time
testSimNPF = cellSimulation(npfGS,npfIC, 10000)
x=collect(1:10000)
plot(x, log1p.(transpose(testSimNPF.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
 savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF2_SplicedCounts.pdf")
 plot(x, log1p.(transpose(testSimNPF.ProteinExpression)), title = "Protein",
  legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
  savefig("C:/Users/rossi/Pictures/CellSimulations/testNPF2_ProteinExpression.pdf")
