include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)

function getNewICS(firstICS, results)

    EpigeneticState = results.EpigeneticState[:,end]
    UnsplicedRNA = results.UnsplicedRNA[:,end]
    SplicingStages = results.SplicingStages
    SplicedRNA = results.SplicedRNA[:,end]
    ProteinExpression = results.ProteinExpression[:,end]

    return (EpigeneticState = EpigeneticState, UnsplicedRNA = UnsplicedRNA,
        SplicingStages= SplicingStages, SplicedRNA = SplicedRNA,
        ProteinExpression = ProteinExpression)
end

function tupleMatrixHcat(tupVec, fieldName::Symbol)
    out = getproperty(tupVec[1], fieldName)
    for i=2:length(tupVec)
        out = hcat(out, getproperty(tupVec[i], fieldName))
    end

    return out
end

function simWithPerturbation(geneSet::Vector{Gene}, ics::NamedTuple, tStepsPer::Int,
     nPerturb::Int; levelToPerturb::String = "ProteinExpression",
     levelToPerturb2::String = "", levelToPerturb3::String = "",
     geneNToPerturb::Int = 1, perturbValue::Int = 0)

    resultsVec = Vector{NamedTuple}(undef,(nPerturb+1))
    for i=1:(nPerturb+1)
        resultsVec[i] = cellSimulation(geneSet, ics, tStepsPer)
        ics = getNewICS(ics, resultsVec[i])
        getproperty(ics, Symbol(levelToPerturb))[geneNToPerturb] = perturbValue
        if levelToPerturb2 == ""
        else
            getproperty(ics, Symbol(levelToPerturb2))[geneNToPerturb] = perturbValue
        end
        if levelToPerturb3 == ""
        else
            getproperty(ics, Symbol(levelToPerturb3))[geneNToPerturb] = perturbValue
        end
    end

    rVec2 = Vector{Matrix{Int}}(undef, length(keys(resultsVec[1])))
    j=1
    for key in keys(resultsVec[1])
        rVec2[j] = tupleMatrixHcat(resultsVec, key)
        j+=1
    end

    return (EpigeneticState = rVec2[1], UnsplicedRNA = rVec2[2],
        SplicingStages= rVec2[3], SplicedRNA = rVec2[4],
        ProteinExpression = rVec2[5])
end

test = simWithPerturbation(testGS, testIC, 100, 9)

#plot results
using Plots
x= collect(1:1000)
plot(x, transpose(test.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_ProteinPert.pdf")

plot(x, log1p.(transpose(test.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_ProteinPert.pdf")

plot(x, log1p.(transpose(test.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_ProteinPert.pdf")

plot(x, log1p.(transpose(test.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_ProteinPert.pdf")

#perturb RNA and protein
sim2 = simWithPerturbation(testGS, testIC, 100, 9, levelToPerturb2 = "SplicedRNA")

plot(x, transpose(sim2.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_ProteinRNAPert.pdf")

plot(x, log1p.(transpose(sim2.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_ProteinRNAPert.pdf")

plot(x, log1p.(transpose(sim2.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_ProteinRNAPert.pdf")

plot(x, log1p.(transpose(sim2.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_ProteinRNAPert.pdf")


#also perturb epigenetic state
sim3 = simWithPerturbation(testGS, testIC, 100, 9, levelToPerturb2 = "SplicedRNA",
    levelToPerturb3 = "EpigeneticState")

plot(x, transpose(sim3.EpigeneticState), title = "Epigenetic State",
 legend = false, xlabel = "Time (minutes)", ylabel = "Accessibile Loci")
savefig("C:/Users/rossi/Pictures/CellSimulations/EpigeneticState_ProteinRNAEpiPert.pdf")

plot(x, log1p.(transpose(sim3.SplicedRNA)), title = "Spliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/SplicedCounts_ProteinRNAEpiPert.pdf")

plot(x, log1p.(transpose(sim3.UnsplicedRNA)), title = "Unspliced Counts",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Counts")
savefig("C:/Users/rossi/Pictures/CellSimulations/UnsplicedCounts_ProteinRNAEpiPert.pdf")

plot(x, log1p.(transpose(sim3.ProteinExpression)), title = "Protein",
 legend = false, xlabel = "Time (minutes)", ylabel = "log1p Molecules")
savefig("C:/Users/rossi/Pictures/CellSimulations/ProteinExpression_ProteinRNAEpiPert.pdf")
