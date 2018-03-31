# illustrative example of using cellAggregator
# after sourcing the functions in cellAggregatorFunctions.R
# last updated 31 March 2018

library(igraph)

timesteps = 100
expression = c(10,2,2,10)
affVal = 3
init = initialise(expression,affVal)
cellInfo = init$cellInfo
proteinInfo = init$proteinInfo
proteins2cell = init$proteins2cell
res = runSimulation(cellInfo,proteinInfo,proteins2cell,
                    timesteps,affVal,
                    expression)

g = graph.adjacency(res$listOfCellNetworks[[1]],mode="undirected")
lyout = layout.circle(g)

par(mfrow=c(5,2))
par(mar=c(0,0,2,0))
for (i in c(1,seq(10,100,length.out=4))){
  print(i)
  plotCellNetwork(res$listOfCellNetworks[[i]],res$cellInfo,layout=lyout)
  title(i)
  plotCellNetwork(res$listOfCellNetworks[[i]],res$cellInfo)
  title(i)
}


plotCellNetwork(res$listOfCellNetworks[[1]],layout=lyout)
for (i in 2:length(res$listOfCellNetworks)) {
plotCellNetwork(res$listOfCellNetworks[[i]],layout=lyout,
                alpha = 0.02,add = TRUE)
}


netedgecols = do.call(rbind,lapply(res$listOfCellNetworks, plotCellNetwork,returnedgecols=TRUE))
par(mfrow=c(2,1))
barplot(t(netedgecols), col = c(colnames(netedgecols)),border = FALSE)
barplot(t(netedgecols/rowSums(netedgecols)), col = c(colnames(netedgecols)),border = FALSE)
