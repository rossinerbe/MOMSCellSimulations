include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testSim = cellSimulation(testGS, testIC, 1000)

testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicedRNA",
 perturbValue = 0, geneNToPerturb = 1, startPerturb = 50,
 levelToPerturb2 = "UnsplicedRNA", perturbValue2 = 0,
  geneNToPerturb2 = 1, startPerturb2 = 50)


using Plots
x= collect(1:1000)
plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t50_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t50_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t50_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t50_ProteinExpression.pdf")


#same for gene7 (highest expressed protein)
testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicedRNA",
 perturbValue = 0, geneNToPerturb = 7, startPerturb = 50,
 levelToPerturb2 = "UnsplicedRNA", perturbValue2 = 0,
  geneNToPerturb2 = 7, startPerturb2 = 50)


plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene7KO_t50_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene7KO_t50_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene7KO_t50_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene7KO_t50_ProteinExpression.pdf")

#KO gene4 (median protein expression)
testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicedRNA",
 perturbValue = 0, geneNToPerturb = 4, startPerturb = 50,
 levelToPerturb2 = "UnsplicedRNA", perturbValue2 = 0,
  geneNToPerturb2 = 4, startPerturb2 = 50)


plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene4KO_t50_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene4KO_t50_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene4KO_t50_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene4KO_t50_ProteinExpression.pdf")


#KO gene 1 after 5tps
testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicedRNA",
 perturbValue = 0, geneNToPerturb = 1, startPerturb = 5,
 levelToPerturb2 = "UnsplicedRNA", perturbValue2 = 0,
  geneNToPerturb2 = 1, startPerturb2 = 5)

plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t5_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t5_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t5_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t5_ProteinExpression.pdf")

#KO gene 1 after 500tps
testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicedRNA",
 perturbValue = 0, geneNToPerturb = 1, startPerturb = 500,
 levelToPerturb2 = "UnsplicedRNA", perturbValue2 = 0,
  geneNToPerturb2 = 1, startPerturb2 = 500)

plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t500_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t500_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t500_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t500_ProteinExpression.pdf")


#KO difference with mRNA degradation allowed, not just instantly 0
testSimKO = cellSimulation(testGS, testIC, 1000, levelToPerturb = "SplicingStages",
 perturbValue = 0, geneNToPerturb = 1, startPerturb = 100)

plot(x, transpose(testSimKO.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t100_wDeg_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimKO.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t100_wDeg_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t100_wDeg_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimKO.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1KO_t100_wDeg_ProteinExpression.pdf")


#make gene 1 constitutively epigenetic repressed (currently it is always accessible)
testSimEpi = cellSimulation(testGS, testIC, 1000, levelToPerturb = "EpigeneticState",
 perturbValue = 0, geneNToPerturb = 1, startPerturb = 0)

plot(x, transpose(testSimEpi.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1EpiRepress_t1_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimEpi.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1EpiRepress_t1_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimEpi.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1EpiRepress_t1_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimEpi.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Gene1EpiRepress_t1_ProteinExpression.pdf")


#degrade all of protein 7 and prevent more being made
testSimProt = cellSimulation(testGS, testIC, 1000, levelToPerturb = "ProteinExpression",
 perturbValue = 0, geneNToPerturb = 7, startPerturb = 50)

plot(x, transpose(testSimProt.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/Protein7_t50_EpigeneticState.pdf")

plot(x, log1p.(transpose(testSimProt.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Protein7_t50_SplicedCounts.pdf")

plot(x, log1p.(transpose(testSimProt.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/Protein7_t50_UnsplicedCounts.pdf")

plot(x, log1p.(transpose(testSimProt.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/Protein7_t50_ProteinExpression.pdf")
