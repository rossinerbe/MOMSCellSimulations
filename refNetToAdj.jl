function edgeListtoAdj(edgeList, orderedNames)
    adjMat = zeros(length(orderedNames), length(orderedNames))
    for i=1:size(edgeList)[1]
        row= findall(x->x == edgeList.Gene1[i], orderedNames)[1]
        col = findall(x->x == edgeList.Gene2[i], orderedNames)[1]
        adjMat[row, col] = 1
    end
     return adjMat
 end
