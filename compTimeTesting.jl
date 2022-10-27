include("simulateCells.jl");

testGS = generateGeneSet(10, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
@time cellSimulation(testGS, testIC, 1000) #0.24s, 1.02M allocation

testGS = generateGeneSet(100, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
@time cellSimulation(testGS, testIC, 1000) #3.83s, 9.18M allocations

testGS = generateGeneSet(1000, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
@time cellSimulation(testGS, testIC, 1000) #280.5s, 88.8M allocations
