using CSV, DataFrames
using Statistics

#get directory paths
directories = Vector{String}(undef,0)
for (root, dirs, files) in walkdir("SimDataRNAVelo")
    for dir in dirs
        push!(directories, joinpath(root, dir)) # path to directories
    end
end

#remove directories to simulation folders
directories = directories[1:12]

#create output data matrix
using TransferEntropy
using Statistics
using StatsBase
using LinearAlgebra
using Plots
include("refNetToAdj.jl");
outRes = Vector{Vector{Int64}}(undef, length(directories))

for k=1:length(directories)
    #read in the data
    exprData = CSV.read(joinpath(directories[k], "ExpressionData.csv"), DataFrame)
    refNet = CSV.read(joinpath(directories[k], "refNetwork.csv"), DataFrame)
    exprDataOnly = exprData[:,2:end]
    exprDataOnly = Matrix(exprDataOnly)

    est =VisitationFrequency(RectangularBinning(200))
    miResTmp = Matrix{Float64}(undef, size(exprData)[1], size(exprData)[1])
    for i=1:size(exprData)[1]
        for j=1:size(exprData)[1]
                miResTmp[i,j] = mutualinfo(exprDataOnly[i,:], exprDataOnly[j,:], est)
        end
    end

    adjMat = edgeListtoAdj(refNet, exprData[:,1])

    rankVec = Vector{Vector{Int}}(undef, size(miResTmp)[2])
    for i=1:size(miResTmp)[2]
        causalGenes = findall(x->x != 0, adjMat[:,i])
        #set diagonal of mi matrix to zero
        diag(miResTmp) .= 0
        ranks = findall(in(causalGenes), reverse(sortperm(miResTmp[:,i])))
        rankVec[i] = ranks
    end

    ranks = vcat(rankVec...)
    rankCounts = countmap(ranks)
    rankCountsOrdered = Vector{Int}(undef, size(miResTmp)[2])
    for i=1:size(miResTmp)[2]
        ind = findall(x->x == i, collect(keys(rankCounts)))
        if length(ind) == 0
            rankCountsOrdered[i] = 0
        else
            rankCountsOrdered[i] = collect(values(rankCounts))[ind[1]]
        end
    end

    outRes[k] = rankCountsOrdered
end

using JLD2
save_object("MIrankCountsBoolODE.jld2", outRes)

bar(outRes[1], orientation=:h, yticks=(1:length(outRes[1])), yflip=true,
 xlabel = "Number of Genes", ylabel = "Mutual information rank", legend = false)
savefig("BoolODE_MIRanks_GSD.pdf")
bar(outRes[12], orientation=:h, yticks=(1:length(outRes[12])), yflip=true,
 xlabel = "Number of Genes", ylabel = "Mutual information rank", legend = false)
savefig("BoolODE_MIRanks_trifurcating.pdf")
