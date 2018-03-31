# last updated 31 March 2018
# cellAggregatorFunctions.R

library(gtools)
library(igraph)
library(parallel)
library(ggplot2)
require(scales)

#####################################################################
#####################################################################
# functions
#####################################################################
# helper functions

initialiseProteinNetwork = function(proteinInfo) {
  proteinNetwork <- matrix(0,nrow=nrow(proteinInfo),ncol=nrow(proteinInfo))
  rownames(proteinNetwork) <- rownames(proteinInfo)
  colnames(proteinNetwork) <- rownames(proteinInfo)
  return(proteinNetwork)
}

cellNetworkFromProteinNetwork <- function(proteinNetwork,proteins2cell) {
  # if any cell in each submatrix is non-zero then set the cellnetwork as 1
  
  cells = names(proteins2cell)
  
  cellNetwork <- matrix(0,nrow=length(proteins2cell),ncol=length(proteins2cell))
  rownames(cellNetwork) <- names(proteins2cell)
  colnames(cellNetwork) <- names(proteins2cell)
  
  for (cell1 in cells) {
    for (cell2 in cells) {
      subMatrix = proteinNetwork[proteins2cell[[cell1]],proteins2cell[[cell2]]]
      if (sum(subMatrix)>0) {
        cellNetwork[cell1,cell2] <- 1
      }
    }
  }
  return(cellNetwork)
}


getproteinStatus = function(proteinNetwork) {
  return(ifelse(rowSums(proteinNetwork)>0,"Bound","Unbound"))
}

speedDateCells = function(cellNetwork) {
  # input is cellNetwork taken from cellNetworkFromProteinNetwork()
  # output is a pairing of cells
  g = graph.adjacency(cellNetwork,mode="undirected")
  positions = layout_with_graphopt(g)
  rownames(positions) <- rownames(cellNetwork)
  d = as.matrix(dist(positions))
  diag(d) <- NA
  speedDate = matrix("",nrow=nrow(cellNetwork)/2,ncol=2)
  for (i in 1:(nrow(speedDate)-1)) {
    
    # randomly select from rownames(d)
    cell1 = sample(rownames(d),1)
    
    # rank of distance, rank 1 is closest
    r = rank(d[cell1,!colnames(d) %in% c(cell1)])
    
    # do 1/distance based probability
    cell2 = sample(names(r),1,prob=(length(r)-r))
    
    speedDate[i,] <- c(cell1,cell2)
    d = d[!rownames(d) %in% c(cell1,cell2),!colnames(d) %in% c(cell1,cell2)]
  }
  speedDate[nrow(speedDate),] <- rownames(d)
  return(speedDate)
}


proteinPairing = function(proteinNetwork,proteinStatus,proteinInfo,proteins2cell,speedDate,affVal) {
  # perform the protein pairing with equal prob to the available proteins
  # if proteins are found to be available then they bind for a 
  # set length of time depending on the protein and its affinity value
  # proteins2cell needs to be defined in the global environment
  # do not allow all speedDating cells to bind at once
  # can lead to weird oscillatory behaviour
  for (i in 1:ceiling(propSpeedDating*nrow(speedDate))) {
    cells = speedDate[i,]
    availableProteins1 = intersect(names(which(proteinStatus=="Unbound")), 
                                   proteins2cell[[cells[1]]])
    availableProteins2 = intersect(names(which(proteinStatus=="Unbound")), 
                                   proteins2cell[[cells[2]]])
    if (length(availableProteins1)==0|length(availableProteins2)==0) next
    protein1 = sample(availableProteins1,1)
    protein2 = sample(availableProteins2,1)
    type1 = as.character(proteinInfo[protein1,1])
    type2 = as.character(proteinInfo[protein2,1])
    if (type1!=type2) {
      proteinNetwork[protein1,protein2] <- diffLength
      proteinNetwork[protein2,protein1] <- diffLength
    }
    if (type1==type2) {
      if (type1=="A") {
        proteinNetwork[protein1,protein2] <- sameLength*affVal
        proteinNetwork[protein2,protein1] <- sameLength*affVal
      } else {
        proteinNetwork[protein1,protein2] <- sameLength
        proteinNetwork[protein2,protein1] <- sameLength
      }
    }
  }
  return(proteinNetwork)
}

countDown = function(proteinNetwork) {
  proteinNetwork = apply(proteinNetwork,1:2,function(x) {
    if (x<=0) return(0) else 
      return(x-1)
  })
  return(proteinNetwork)
}

mixingIndex = function(cellNetwork,cellInfo) {
 
  cellNetwork2 = cellNetwork
    
  prop = NULL
  for (cell in rownames(cellNetwork2)) {
    cellcol = as.character(cellInfo[cell,1])
    nms = names(which(cellNetwork2[cell,]==1))
    if (length(nms)==0) {
      prop[cell] <- NA
    } else {
    prop[cell] <- sum(as.character(cellInfo[nms,1])!=cellcol)/length(nms)
      }
  }
  
  greengreen = sum(cellNetwork[cellInfo$Colour=="green",
                                   cellInfo$Colour=="green"])
  greenred = sum(cellNetwork[cellInfo$Colour=="green",
                                 cellInfo$Colour=="red"])
  redgreen = sum(cellNetwork[cellInfo$Colour=="red",
                                 cellInfo$Colour=="green"])
  redred = sum(cellNetwork[cellInfo$Colour=="red",
                               cellInfo$Colour=="red"])
  N = greengreen+greenred+redgreen+redred
  mixN = greenred + redgreen
  if (N==0) {
    edgeMixingIndex <- NA
  } else {
    edgeMixingIndex = mixN/N
  }
  return(list(meanMixingIndex = mean(prop,na.rm=TRUE),
              sdMixingIndex = sd(prop,na.rm=TRUE),
              edgeMixingIndex = edgeMixingIndex,
              prop=prop))
}

#####################################################################
# functions to initialise

initialise <- function(expression, affVal) {
  
  cellColour <- c(rep("green",numCellsGreen),rep("red",numCellsRed))
  cellProtA <- c(rep(expression[1],numCellsGreen),rep(expression[3],numCellsRed))
  cellProtB <- c(rep(expression[2],numCellsGreen),rep(expression[4],numCellsRed))
  
  cellInfo <- data.frame(
    Colour = cellColour,
    NumberProteinA = cellProtA,
    NumberProteinB = cellProtB
  )
  rownames(cellInfo) <- paste0("cell_",sprintf("%03d",1:nrow(cellInfo)))
  
  proteinType <- NULL
  proteinCell <- NULL
  proteinCellColour <- NULL
  
  for (cell in rownames(cellInfo)) {
    numProtein <- sum(cellInfo[cell,2:3])
    proteinType <- c(proteinType,rep(c("A","B"),cellInfo[cell,2:3]))
    proteinCell <- c(proteinCell,rep(cell,numProtein))
    proteinCellColour <- c(proteinCellColour,rep(as.character(cellInfo[cell,1]),numProtein))
  }
  
  proteinInfo <- data.frame(
    proteinType <- proteinType,
    proteinCell <- proteinCell,
    proteinCellColour <- proteinCellColour
  )
  rownames(proteinInfo) <- paste0("protein_",sprintf("%03d",1:nrow(proteinInfo)))
  
  proteins2cell <- split(rownames(proteinInfo),proteinInfo[,2])
  
  return(list(cellInfo=cellInfo,
              proteinInfo=proteinInfo,
              proteins2cell=proteins2cell))
}

#####################################################################
# function to run simulation
# usage:
# expression = c(3,3,3,3) # green A and B, red A and B
# affVal = 1
# initialise(expression,affVal)

runSimulation <- function(cellInfo,proteinInfo,proteins2cell,
                          timesteps,affVal,expression,
                          verbose=FALSE) {
  
  # input
  # timesteps
  # affVal
  # expression
  # proteinInfo
  # proteins2cell
  # cellInfo
  
  propMixingIndex = matrix(NA,ncol=timesteps,nrow=nrow(cellInfo))
  meanMixingIndex = NULL
  sdMixingIndex = NULL
  edgeMixingIndex = NULL
  listOfCellNetworks = list()
  
  for (timestep in 1:timesteps) {
    if (verbose) print(timestep)
    if (timestep==1) {
      proteinNetwork = initialiseProteinNetwork(proteinInfo)
    }
    proteinStatus = getproteinStatus(proteinNetwork)
    cellNetwork = cellNetworkFromProteinNetwork(proteinNetwork,proteins2cell)
    listOfCellNetworks[[timestep]] <- cellNetwork
    mi = mixingIndex(cellNetwork,cellInfo)
    meanMixingIndex[timestep] <- mi$meanMixingIndex
    sdMixingIndex[timestep] <- mi$sdMixingIndex
    edgeMixingIndex[timestep] <- mi$edgeMixingIndex
    if (!is.null(mi$prop)) propMixingIndex[,timestep] <- mi$prop
    if (timestep==timesteps) break
    speedDate = speedDateCells(cellNetwork)
    proteinNetwork <- proteinPairing(proteinNetwork,proteinStatus,proteinInfo,proteins2cell,speedDate,affVal)
    proteinNetwork <- countDown(proteinNetwork)
  }
  
  # unchanged/hardcoded parameters:
  diffLength <- 1 # length of timesteps that different proteins stay bound
  sameLength <- 3 # length of timesteps that same proteins stay bound
  numCellsGreen <- 25
  numCellsRed <- 25
  timesteps <- 100 # number of timesteps in simulation
  propSpeedDating <- 0.75 # proportion of speed dating cell pairs to allow to bind
  propTimestepsAverage <- 0.75 # proportion of timesteps i average over to get the statistics (1-this value is the 'burn in' proportion)
  
  inputParameters = c("diffLength"=diffLength,
                      "sameLength"=sameLength,
                      "numCellsGreen"=numCellsGreen,
                      "numCellsRed"=numCellsRed,
                      "timesteps"=timesteps,
                      "propSpeedDating"=propSpeedDating,
                      "propTimestepsAverage"=propTimestepsAverage)
  
  # values returned are the last timesteps cellNetwork
  # the meanMixingIndex averaged over the last half of timesteps
  # the sdMixingIndex averaged over the last half of timesteps
  
  return(list(
    inputParameters = inputParameters,
    affVal = affVal,
    expression = expression,
    cellInfo = cellInfo,
    cellNetwork=cellNetwork,
    mean = mean(meanMixingIndex[(ceiling((1-propTimestepsAverage)*timesteps):timesteps)]),
    sd = mean(sdMixingIndex[(ceiling((1-propTimestepsAverage)*timesteps):timesteps)]),
    edge = mean(edgeMixingIndex[(ceiling((1-propTimestepsAverage)*timesteps):timesteps)]),
    listOfCellNetworks = listOfCellNetworks,
    propMixingIndex=propMixingIndex
  ))
}

simulationWrap = function(expression,affVal,...) {
  init = initialise(expression,affVal)
  cellInfo = init$cellInfo
  proteinInfo = init$proteinInfo
  proteins2cell = init$proteins2cell
  res = runSimulation(cellInfo,proteinInfo,proteins2cell,
                      timesteps,affVal,
                      expression,...)
  return(res)
}

#####################################################################
# plotting functions

plotNetworkIgraph = function(res) {
  cellNetwork = res$cellNetwork
  cellInfo = res$cellInfo
  g = graph.adjacency(cellNetwork,mode="undirected")
  V(g)$color = as.character(cellInfo$Colour)
  V(g)$size = 7
  plot(g, layout=layout_with_graphopt)
}

plotCellNetwork = function(cellNetwork,cellInfo,labels=FALSE,layout=NULL,points=FALSE,
                           alpha = 1,returnedgecols=FALSE,...) {
  g = graph.adjacency(cellNetwork,mode="undirected")
  V(g)$color = as.character(cellInfo$Colour)
  cols = V(g)$color
  names(cols) <- names(V(g))
  edgecols = apply(t(apply(get.edgelist(g),1,function(x)cols[x])),1,function(x) {
    if (length(unique(x))==1) {
      return(x[1])
    } else {
        return("orange")
      }
  })
  E(g)$color = edgecols
  if (returnedgecols) return(table(factor(edgecols,levels=c("green","red","orange"))))
  if (!labels) {
    V(g)$label = ""
  }
  V(g)$size = 7
  if (is.null(layout)) layout = layout_with_graphopt(g)
  plot(g, layout=layout, edge.color = alpha(edgecols, alpha),...)
}

convertSummaryValuesToExpressionVector = function(summaryExpressionRatio,
                                                  summaryExpressionScenario) {
  # examples
  # summaryExpressionScenario = "Ab_B"
  # summaryExpressionRatio = 5
  summaryExpressionRatio = as.numeric(summaryExpressionRatio)
  
  summaryExpressionScenarioSplit = unlist(strsplit(summaryExpressionScenario,""))
  expression = rep(0,4)
  expression[summaryExpressionScenarioSplit %in% c("A","B")] = summaryExpressionRatio
  expression[summaryExpressionScenarioSplit %in% c("a","b")] = 1
  
  return(expression)
}
