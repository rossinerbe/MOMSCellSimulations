include("simulateCells.jl");

simData = Array{Float64}(undef, 10, 2000, 1000)
testGS = generateGeneSet(10, allowFeedbackLoops = true)
for i=1:2000
    testIC = generateInitialConditions(testGS, seed = i)
    testSim = cellSimulation(testGS, testIC, 1000, seed = i)
    #get total counts
    simData[:,i,:] = testSim.UnsplicedRNA + testSim.SplicedRNA
end

t1 = Float32.(log1p.(simData[:,:,101]))
t2 = Float32.(log1p.(simData[:,:,102]))

#train RNAForecaster
include("trainRNAForecaster.jl")
include("makeRecursivePredictions.jl")

trainingResults = createEnsembleForecaster(t1, t2,
 nNetworks = 10, checkStability = false, hiddenLayerNodes = 100)

#now see what the loss is over simulated time as we try to recursively
#predict expression into the future states on which the model was not
#trained on
exprPreds = ensembleExpressionPredictions(trainingResults, t1, 50)


#compare to simulated benchmark
predictionErrors = Vector{Float32}(undef, size(exprPreds)[3])
for j=1:size(exprPreds)[3]
    actual = Float32.(log1p.(simData[:,:,(j+101)]))
    #calculate average cell-wise mse
    predictionErrors[j] = mse(exprPreds[:,:,j], actual, agg= sum)/(size(actual)[1] * size(actual)[2])

end

using JLD2
save_object("RNAForecasterPredictionErrors1.jld2", predictionErrors)

predictionErrors = load_object("RNAForecasterPredictionErrors1.jld2")
using Plots
x=collect(1:50)
plot(x,log1p.(predictionErrors), legend = false)
savefig("C:/Users/rossi/Pictures/CellSimulations/RNAForecasterlogMSE_50tps.pdf")
