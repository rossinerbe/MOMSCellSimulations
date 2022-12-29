# Multi-Omic Mechanistic Simulations (MOMS) 
Example to run a random 100 gene simulation for 1000 time points for one gene:
 ```
include("simulateCells.jl");

testGS = generateGeneSet(100, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)

testSim = cellSimulation(testGS, testIC, 1000)
```

To create a data set with the same number of genes spanning 1000 cells we instead call:
```
testGS = generateGeneSet(100, allowFeedbackLoops = true)
testIC = generateInitialConditions(testGS)
testMCSim = makeMultiCellDataSet(testGS, testIC, 1000, 1000)
```
To get out the adjacency matrix of TF-gene realtionships that paramterize the simulation we can use the function:
```
adjMat = getTFGeneAMat(testGS)
```
