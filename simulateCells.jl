#unit time in this simulation is going to be based on time to transcribe/translate a gene,
#which are helpfully very similar. Assumed to be ~1 min for our purposes

struct Gene
   Name::String
   #rate/probability of transcription when accessible
   BaseTranscriptionProbability::Float64
   #which genes increase transcription if accessible
   TranscriptionActivators::Vector{String}
   #which genes inhibit transcription
   TranscriptionRepressors::Vector{String}
   #which genes can switch gene state from repressed accessible
   EpigeneticUpregulators::Vector{String}
    #which genes can switch gene state from accessible to repressed
   EpigeneticDownregulators::Vector{String}
   SplicingRate::Int
   RNAHalfLife::Int
   ProteinHalfLife::Int
   #probability of unbound mRNA being bound by ribosome and translated/minute
   TranslationInitiationProbability::Float64
   #which genes can bind RNA and prevent translation
   TranslationalInhibitor::Vector{String}
   #which genes can bind protein and cause it to be degraded
   ProteinDegradationFactors::Vector{String}
end


using Random, StatsBase
using Distributions
#function to generate a random set of genes
function generateGeneSet(nGenes::Int; BTProbRange::Vector{Float64} = [0.01,0.001], #proability in a given minute
  maxTAs::Int = 2, maxTRs::Int = 1,
  maxERU::Int = 2, maxERD::Int = 2,
  splicingRateRange::Vector{Int} = [5,10], #http://book.bionumbers.org/what-is-faster-transcription-or-translation/
  maxTI::Int = 1,
  RNAHLRange::Vector{Int} = [60, 900], #http://book.bionumbers.org/how-fast-do-rnas-and-proteins-degrade/
  ProteinHLRange::Vector{Int} = [720, 3600],
  TIPRange::Vector{Float64} = [0.1, 0.75], #https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007070
  maxPDFs::Int = 2,
  allowFeedbackLoops::Bool = false,
  geneNames::Vector{String} = Vector{String}(undef,0),
  seed::Int = 123)

  if length(geneNames) == 0
      geneNames = Vector{String}(undef, nGenes)
      for i=1:nGenes
          geneNames[i] = "Gene" * string(i)
      end
  end

  #generate gene attributes
  geneVec = Vector{Gene}(undef,nGenes)
  Random.seed!(seed)
  for i=1:nGenes
      tmpName = geneNames[i]
      if allowFeedbackLoops
              otherGeneNames = geneNames
      else
              otherGeneNames = geneNames[setdiff(1:end,i)]
      end

              #assign base transcription rate probability per unit time (~1min)
      btrProb = rand(Uniform(BTProbRange[2], BTProbRange[1]))

      #which genes' products are transcriptional activators/repressors of this gene?
      nTAs = rand(1:maxTAs)
      nTRs = rand(0:maxTRs)
      tmpTregs = geneNames[sample(1:length(geneNames),(nTAs+nTRs), replace = false)]
      tmpTAs = tmpTregs[1:nTAs]
      if nTRs > 0
              tmpTRs = tmpTregs[(nTAs+1):end]
      else
              tmpTRs = Vector{String}(undef,0)
      end

      #which genes' products epigenetically upregulate/downregulate it
      nERU = rand(1:maxERU)
      nERD = rand(1:maxERD)
      tmpEregs = geneNames[sample(1:length(geneNames),(nERU+nERD), replace = false)]
      tmpERUs = tmpEregs[1:nERU]
      tmpERDs = tmpEregs[(nERU+1):end]

      #splicing rate
      tmpSR = rand(splicingRateRange[1]:splicingRateRange[2])

      #RNA half life
      tmpRNAHL = rand(RNAHLRange[1]:RNAHLRange[2])

      #protein half life
      tmpProteinHL = rand(ProteinHLRange[1]:ProteinHLRange[2])

      #translation initiation probability
      tmpTIP = rand(Uniform(TIPRange[1], TIPRange[2]))

      #which genes' products inhibit its translation?
      nTI = rand(0:maxTI)
      tmpTIs = geneNames[sample(1:length(geneNames),nTI, replace = false)]

      #which are protein degradation factors
      nPDFs = rand(1:maxPDFs)
      tmpPDFs = otherGeneNames[sample(1:length(otherGeneNames),nPDFs, replace = false)]

      geneVec[i] = Gene(tmpName, btrProb, tmpTAs, tmpTRs, tmpERUs, tmpERDs,
                tmpSR, tmpRNAHL, tmpProteinHL, tmpTIP, tmpTIs, tmpPDFs)
  end

  return(geneVec)
end


function bindingProb(concentration)
        1 .- exp.(-1 .* concentration)
end

#expected degradation in one unit time based on half life formula
function expectedDegradation(halfLife, initialQuantity)
        initialQuantity - (initialQuantity * (0.5^(1/halfLife)))
end

#needs epignetic state, spliced and unspliced counts (with splicing progress),
# and protein expression
function generateInitialConditions(geneSet::Vector{Gene}; maxUnsplicedRNA::Int = 2,
        maxSplicedRNA::Int = 2, maxProtein::Int = 2, seed::Int = 123)

        Random.seed!(seed)

        #get max splicing time for splicing progress tracking matrix
        splicingRates = Vector{Int}(undef, length(geneSet))
        for i=1:length(geneSet)
                splicingRates[i] = geneSet[i].SplicingRate
        end
        maxSR = maximum(splicingRates)

        ics = (EpigeneticState = Vector{Int}(undef, length(geneSet)),
               UnsplicedRNA = Vector{Int}(undef, length(geneSet)),
               SplicingStages = Matrix{Int}(undef, length(geneSet), maxSR+1),
               SplicedRNA = Vector{Int}(undef, length(geneSet)),
               ProteinExpression = Vector{Int}(undef, length(geneSet)))

         ics.SplicingStages .= 0

        for i=1:length(geneSet)
                ics.EpigeneticState[i] = rand(0:2)
                ics.UnsplicedRNA[i] = rand(0:maxUnsplicedRNA)
                ics.SplicingStages[i,end] = ics.UnsplicedRNA[i]
                ics.SplicedRNA[i] = rand(0:maxSplicedRNA)
                ics.ProteinExpression[i] = rand(0:maxProtein)
        end

        return(ics)

end

#simulate a cell with these gene relationships
#might need to eventually add in a variable for total number of ribosomes to account for the
#fact that they can become saturated and slow translation


#the impact an additional TF molecule has on the probability that a corresponding locus is bound by said TF
#should be proportional to the concentration of the TF in the nucleus
#thus, the volume of the nucleus is import. Normal nuclear diameter ranges between
#2 and 10 microns, according to http://book.bionumbers.org/how-big-are-nuclei/
#which yields a volume range of ~4 - 525 um^3

function cellSimulation(genes::Vector{Gene}, initialConditions, timePoints::Int;
    epigeneticLeakiness::Float64 = 0.05, nuclearVolume::Float64 = 125.0,
    cellVolume::Float64 = 2000.0, seed::Int = 123,
    levelToPerturb::String = "None", perturbValue::Int = 0, geneNToPerturb::Int = 1,
    startPerturb::Int = 50, levelToPerturb2::String = "None", perturbValue2::Int = 0,
    geneNToPerturb2::Int = 2, startPerturb2::Int = 50)

    Random.seed!(seed)
    cytosolVolume = cellVolume - nuclearVolume
    geneNames = Vector{String}(undef, length(genes))
    for i=1:length(genes)
            geneNames[i] = genes[i].Name
    end

    ##create output results object
    #contains matrix of epigenetic states over time
    #matrix of unspliced RNA over time
    #matrix tracking splicing progress of immature RNA molecules
    #matrix of spliced RNA over time
    #matrix of protein expression

    #get max splicing time for splicing progress tracking matrix
    splicingRates = Vector{Int}(undef, length(genes))
    for i=1:length(genes)
            splicingRates[i] = genes[i].SplicingRate
    end
    maxSR = maximum(splicingRates)

    results = (EpigeneticState = Matrix{Int}(undef, length(genes), timePoints),
                UnsplicedRNA = Matrix{Int}(undef, length(genes), timePoints),
                SplicingStages = Matrix{Int}(undef, length(genes), maxSR+1),
                SplicedRNA = Matrix{Int}(undef, length(genes), timePoints),
                ProteinExpression = Matrix{Int}(undef, length(genes), timePoints))

    #add initial conditions to results object
    results.EpigeneticState[:,1] = initialConditions.EpigeneticState
    results.UnsplicedRNA[:,1] = initialConditions.UnsplicedRNA
    results.SplicingStages[:,:] = initialConditions.SplicingStages
    results.SplicedRNA[:,1] = initialConditions.SplicedRNA
    results.ProteinExpression[:,1] = initialConditions.ProteinExpression

    for i=1:(timePoints-1)
        for j=1:length(genes)
            ##RNA transcription
            es = results.EpigeneticState[j,i]
            transcriptCountChange = [0,0]
            if es==2
                    lociRepressed = [0,0]
            elseif es == 1
                    lociRepressed = [1,0]
            elseif es ==0
                    lociRepressed = [1,1]
            end
            #calculate new transcripts produced
            #determine whether repressor binds to each locus, then no transcription
            if lociRepressed[1] == 1 && lociRepressed[2] == 1
            else
                    if length(findall(in(genes[j].TranscriptionRepressors), geneNames)) > 0
                            repressorConcs = results.ProteinExpression[findall(in(genes[j].TranscriptionRepressors), geneNames), i] ./ nuclearVolume
                            bindingProbsR = bindingProb(repressorConcs)
                            randResults = bindingProbsR .> rand(length(bindingProbsR))
                            if length(findall(x->x==1, randResults)) > 0
                                    lociRepressed[1] = 1
                            end

                            #other genetic loci
                            randResults = bindingProbsR .> rand(length(bindingProbsR))
                            if length(findall(x->x==1, randResults)) > 0
                                    lociRepressed[2] = 1
                            end
                    end
                    #determine whether activator binds to each locus; then transcription occurs
                    if lociRepressed[1] == 1 && lociRepressed[2] == 1
                    elseif lociRepressed[1] == 0
                            if length(findall(in(genes[j].TranscriptionActivators), geneNames)) > 0
                                    activatorConcs = results.ProteinExpression[findall(in(genes[j].TranscriptionActivators), geneNames), i] ./ nuclearVolume
                                    bindingProbsA = bindingProb(activatorConcs)
                                    randResults = bindingProbsA .> rand(length(bindingProbsA))
                                    if length(findall(x->x==1, randResults)) == 1
                                         transcriptCountChange[1] = 1
                                    end
                            end
                    elseif lociRepressed[2] == 0
                            if length(findall(in(genes[j].TranscriptionActivators), geneNames)) > 0
                                    activatorConcs = results.ProteinExpression[findall(in(genes[j].TranscriptionActivators), geneNames), i] ./ nuclearVolume
                                    bindingProbsA = bindingProb(activatorConcs)
                                    randResults = bindingProbsA .> rand(length(bindingProbsA))
                                    if length(findall(x->x==1, randResults)) == 1
                                         transcriptCountChange[2] = 1
                                    end
                            end
                    end
                    #if neither occurs, do basic TFs bind and lead to transcription regardless?
                    if lociRepressed[1] == 0 && transcriptCountChange[1] == 0
                            randResults = genes[j].BaseTranscriptionProbability > rand()
                            if randResults == true
                                 transcriptCountChange[1] = 1
                            end
                    end
                    if lociRepressed[2] == 0 && transcriptCountChange[2] == 0
                            randResults = genes[j].BaseTranscriptionProbability > rand()
                            if randResults == true
                                 transcriptCountChange[2] = 1
                            end
                    end
            end
            #if not accessible does transcription occur anyway (lower probability), leakiness modifier
            if lociRepressed[1] == 1
                    randResults = (genes[j].BaseTranscriptionProbability * epigeneticLeakiness) > rand()
                    if randResults == true
                         transcriptCountChange[1] = 1
                    end
            end
            if lociRepressed[2] == 1
                    randResults = (genes[j].BaseTranscriptionProbability * epigeneticLeakiness) > rand()
                    if randResults == true
                         transcriptCountChange[2] = 1
                    end
            end

            ##RNA degradation
            exDeg = expectedDegradation(genes[j].RNAHalfLife, results.SplicedRNA[j,i])
            #separate into integer and decimal portions
            decPart = exDeg - round(exDeg, RoundDown)
            intPart = Int(round(exDeg, RoundDown))
            #for decimal part, use random number generator to determine if the transcript is degraded
            if decPart > rand()
                    intPart += 1
            end

            ##RNA splicing
            newSpliced = results.SplicingStages[j,end]
            #move all unspliced transcripts one minute further in splicing process
            results.SplicingStages[j,2:end] = results.SplicingStages[j,1:(size(results.SplicingStages)[2]-1)]
            results.SplicingStages[j,1] = 0
            results.SplicingStages[j,(size(results.SplicingStages)[2] - genes[j].SplicingRate)] = sum(transcriptCountChange)

            #set the next time points of RNA expression, spliced and unspliced
            results.UnsplicedRNA[j,i+1] = results.UnsplicedRNA[j,i] + sum(transcriptCountChange) - newSpliced
            results.SplicedRNA[j,i+1] = results.SplicedRNA[j,i] + newSpliced - intPart

            ##Protein Translation
            newProteins = 0
            if results.SplicedRNA[j,i] == 0
                    #if theres no RNA no translation
            elseif sum(results.ProteinExpression[findall(in(genes[j].TranslationalInhibitor), geneNames),i]) > 0
                    inhibitorConcs = results.ProteinExpression[findall(in(genes[j].TranslationalInhibitor), geneNames),i] ./ cytosolVolume
                    bindingProbsI = bindingProb(inhibitorConcs)
                    unboundRNAs = results.SplicedRNA[j,i]
                    for k=1:length(bindingProbsI)
                           randBindingRes = rand(unboundRNAs)
                           unboundRNAs -= length(findall(x->x < bindingProbsI[k], randBindingRes))
                   end

                   probVec = rand(unboundRNAs)
                   newProteins = newProteins + sum(genes[j].TranslationInitiationProbability .> probVec)

            else
                    availableRNACount = results.SplicedRNA[j,i]
                    probVec = rand(availableRNACount)
                    newProteins = newProteins + sum(genes[j].TranslationInitiationProbability .> probVec)
            end

            ##Protein Degradation
            #normal expectedDegradation
            exDeg = expectedDegradation(genes[j].ProteinHalfLife, results.ProteinExpression[j,i])
            #separate into integer and decimal portions
            decPart = exDeg - round(exDeg, RoundDown)
            intPart = Int(round(exDeg, RoundDown))
            #for decimal part, use random number generator to determine if the transcript is degraded
            if decPart > rand()
                    intPart += 1
            end

            #results.ProteinExpression[j,i+1] = results.ProteinExpression[j,i] + newProteins - intPart

            #degradation caused by protein-specific factors
            if results.ProteinExpression[j,i] == 0
                    #if theres no protein, there's no degradation of it either
                    #update protein levels at next time point
                    results.ProteinExpression[j,i+1] = newProteins
            elseif sum(results.ProteinExpression[findall(in(genes[j].ProteinDegradationFactors), geneNames),i]) > 0
                    degConcs = results.ProteinExpression[findall(in(genes[j].ProteinDegradationFactors), geneNames),i] ./ cytosolVolume
                    bindingProbsD = bindingProb(degConcs)
                    unboundProteins = results.ProteinExpression[j,i]
                    for k=1:length(bindingProbsD)
                            if unboundProteins < 0
                                    continue
                            end
                           randBindingRes = rand(unboundProteins)
                           unboundProteins -= length(findall(x->x < bindingProbsD[k], randBindingRes))
                   end
                   degradationEffect = unboundProteins - intPart
                   if degradationEffect < 0
                           degradationEffect = 0
                   end

                   results.ProteinExpression[j,i+1] = degradationEffect + newProteins
           else
                   degradationEffect = results.ProteinExpression[j,i] - intPart
                   if degradationEffect < 0
                           degradationEffect = 0
                   end
                   results.ProteinExpression[j,i+1] = degradationEffect + newProteins
           end

            ##Epigenetic state changes
            euCount = sum(results.ProteinExpression[findall(in(genes[j].EpigeneticUpregulators), geneNames), i])
            edCount = sum(results.ProteinExpression[findall(in(genes[j].EpigeneticDownregulators), geneNames), i])
            epiState = results.EpigeneticState[j,i]
            if  epiState == 0 && euCount > 0
                    #chance of increased accessibility at both loci
                    activatorConcs = results.ProteinExpression[findall(in(genes[j].EpigeneticUpregulators), geneNames), i] ./ nuclearVolume
                    bindingProbsA = bindingProb(activatorConcs)
                    for k=1:2
                            if sum(bindingProbsA .> rand(length(bindingProbsA))) > 1
                                    epiState +=1
                            end
                    end
            elseif epiState == 1 && euCount > 0
                    #chance of increased accessibility at one loci
                    activatorConcs = results.ProteinExpression[findall(in(genes[j].EpigeneticUpregulators), geneNames), i] ./ nuclearVolume
                    bindingProbsA = bindingProb(activatorConcs)
                    if sum(bindingProbsA .> rand(length(bindingProbsA))) > 1
                            epiState +=1
                    end
            elseif epiState == 1 && edCount > 0
                    #chance of decreased accessibility at one loci
                    repressorConcs = results.ProteinExpression[findall(in(genes[j].EpigeneticDownregulators), geneNames), i] ./ nuclearVolume
                    bindingProbsR = bindingProb(repressorConcs)
                    if sum(bindingProbsR .> rand(length(bindingProbsR))) > 1
                            epiState -= 1
                    end
            elseif epiState == 2 && edCount > 0
                    #chance of decreased accessibility at both loci
                    repressorConcs = results.ProteinExpression[findall(in(genes[j].EpigeneticDownregulators), geneNames), i] ./ nuclearVolume
                    bindingProbsR = bindingProb(repressorConcs)
                    for k=1:2
                            if sum(bindingProbsR .> rand(length(bindingProbsR))) > 1
                                    epiState -= 1
                            end
                    end
            end

            #set next epigenetic state
            results.EpigeneticState[j,i+1] = epiState
            end
            #perturbation

            if levelToPerturb == "None"
            elseif levelToPerturb == "SplicingStages"
                    if i >= startPerturb
                            getproperty(results, Symbol(levelToPerturb))[geneNToPerturb,:] .= perturbValue
                    end
            else
                    if i >= startPerturb
                            getproperty(results, Symbol(levelToPerturb))[geneNToPerturb, i+1] = perturbValue
                    end
            end

            if levelToPerturb2 == "None"
            elseif levelToPerturb2 == "SplicingStages"
                    if i >= startPerturb2
                            getproperty(results, Symbol(levelToPerturb2))[geneNToPerturb2,:] .= perturbValue2
                    end
            else
                    if i >= startPerturb2
                            getproperty(results, Symbol(levelToPerturb2))[geneNToPerturb2, i+1] = perturbValue2
                    end
            end
        end

        return(results)
end

#get adjacency matrix of TF-gene relationships
function getTFGeneAMat(geneSet::Vector{Gene})
        nGenes = length(geneSet)
        rnaAMat = Int.(zeros(nGenes, nGenes))
        geneNames = Vector{String}(undef, nGenes)
        for i=1:nGenes
           geneNames[i] = "Gene" * string(i)
        end

        for i=1:nGenes
            rnaAMat[findall(in(geneSet[i].TranscriptionActivators), geneNames),i] .= 1
            rnaAMat[findall(in(geneSet[i].TranscriptionRepressors), geneNames),i] .= -1
        end

        return(rnaAMat)
end


#get adjacency matrices from gene set
function getMultilevelRN(geneSet::Vector{Gene})
    nGenes = length(geneSet)
    epgAMat = Int.(zeros(nGenes, nGenes))
    rnaAMat = Int.(zeros(nGenes, nGenes))
    transInhAMat = Int.(zeros(nGenes, nGenes))
    protDAMat = Int.(zeros(nGenes, nGenes))

    geneNames = Vector{String}(undef, nGenes)
   for i=1:nGenes
       geneNames[i] = "Gene" * string(i)
   end

    for i=1:nGenes
        rnaAMat[findall(in(geneSet[i].TranscriptionActivators), geneNames),i] .= 1
        rnaAMat[findall(in(geneSet[i].TranscriptionRepressors), geneNames),i] .= -1

        epgAMat[findall(in(geneSet[i].EpigeneticUpregulators), geneNames),i] .= 1
        epgAMat[findall(in(geneSet[i].EpigeneticDownregulators), geneNames),i] .= -1

        transInhAMat[findall(in(geneSet[i].TranslationalInhibitor), geneNames),i] .= -1

        protDAMat[findall(in(geneSet[i].ProteinDegradationFactors), geneNames),i] .= -1
    end

    return (Epigenetic = epgAMat, RNA = rnaAMat, TranslationInhibition = transInhAMat,
            ProteinDegradation = protDAMat)
end


function makeMultiCellDataSet(genes::Vector{Gene}, initialConditions,
        timePoints::Int, nCells::Int;
            epigeneticLeakiness::Float64 = 0.05, nuclearVolume::Float64 = 125.0,
            cellVolume::Float64 = 2000.0, seed::Int = 123,
            levelToPerturb::String = "None", perturbValue::Int = 0, geneNToPerturb::Int = 1,
            startPerturb::Int = 50, levelToPerturb2::String = "None", perturbValue2::Int = 0,
            geneNToPerturb2::Int = 2, startPerturb2::Int = 50)

          out = Vector{NamedTuple}(undef, nCells)
          for i=1:nCells
                out[i] = cellSimulation(genes, initialConditions,
                        timePoints, epigeneticLeakiness = epigeneticLeakiness,
                        nuclearVolume = nuclearVolume, cellVolume = cellVolume, seed = seed+i,
                        levelToPerturb = levelToPerturb, perturbValue = perturbValue,
                        geneNToPerturb = geneNToPerturb, startPerturb = startPerturb,
                        levelToPerturb2 = levelToPerturb2, perturbValue2 = perturbValue2,
                        geneNToPerturb2 = geneNToPerturb2, startPerturb2 = startPerturb2)
        end
        return(out)
end
