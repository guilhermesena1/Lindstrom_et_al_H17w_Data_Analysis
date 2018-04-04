## ------------------------------------------------------------------------
## Transpose seurat data into WGCNA object
GetDatExpr <- function(srt, markers) {
  genes.use <- unique(markers$gene)
  t(srt@scale.data[genes.use,])
}

## ------------------------------------------------------------------------
# Picks soft threshold power
getSft <- function(datExpr, plot=F) {
  powers <- c(1:10, 2*6:10)
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
}

## ------------------------------------------------------------------------
# Finds modules from raw expression
Modules <- function(datExpr, p) {
  blockwiseModules(datExpr, power = p, TOMType="signed", 
                   minModuleSize = 10, deepSplit = 3, pamRespectsDendro = F,  
                   numericLabels = F, minKMEtoStay = 0, saveTOMs = F, verbose = 5)
}

## ------------------------------------------------------------------------
# Plots Cluster Dendrogram
PlotNet <- function(net) {
  plotDendroAndColors(net$dendrograms[[1]], net$unmergedColors[net$blocks == 1], dendroLabels = F, "Module colors")
}

## ------------------------------------------------------------------------
# Heatmap of a specific module
PlotModule <- function(srt, datExpr, net, which) {
  genes <- colnames(datExpr)[net$colors == which]
  if(length(genes) > 200){
    genes <- genes[1:200]
  }
  DoHeatmap(srt, genes.use = genes, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE, draw.line = T)
}

## ------------------------------------------------------------------------
## Prints top genes (based on connectivity) for a given module
GetGenesFromModule <- function(datExpr, power, net, which) {
  datExpr <- datExpr[,net$colors == which]
  adj <- adjacency(datExpr, power = power)
  adj <- rowSums(adj)
  
  colnames(datExpr)[order(-adj)]
}

## ------------------------------------------------------------------------
## Draws module network using igraph
# The thresh parameter is a hard TOM value cutoff for two genes
# to be connected
PlotModuleNetwork <- function (datExpr, power, net, which, thresh = .01) {
  dat <- datExpr[, net$colors == which]
  IMConn <- softConnectivity(dat)
  ntop <- min(150, sum(net$colors == which))
  top <- rank(-IMConn) <= ntop

  adj <- adjacency(dat[,top], power = power)
  adj[adj <= thresh] <- 0
  adj[adj >= thresh] <- 1
  diag(adj) <- 0

  gg <- simplify(graph_from_adjacency_matrix(adj, mode="undirected"))
  gg <- delete.vertices(gg, which(degree(gg) == 0)) 
  plot(gg, vertex.size = 4, vertex.label.cex = .8, vertex.label.dist = .3, 
       vertex.color="red", vertex.frame.color="black")
}

## ------------------------------------------------------------------------
## Calculates Module Score for all cells
ModuleScore <- function(datExpr, power, net, mon) {
    ans <- net$MEs
    rownames(ans) <- rownames(datExpr)
    
    ans$Pseudotime <- mon@phenoData@data[rownames(ans),]$Pseudotime
    ans <- ans[apply(ans, 1, function(x) {sum(is.na(x)) == 0}),]
    ans
}

## ------------------------------------------------------------------------
## Plots scores
PlotScores <- function(scores, net) {
  #X11(width = 16, height = 9)
  plot(NA, xlab = "Pseudotime", ylab="Module score", 
       xlim = c(0, max(scores$Pseudotime)), ylim = c(min(scores[, !(colnames(scores) %in% "Pseudotime")]), max(scores[,!(colnames(scores) %in% "Pseudotime")])))
  
  colors <- c()
  for(i in unique(net$colors)) {
    if(i != "grey"){
        if(sum(net$colors == i)>10) {
          lines(smooth.spline(scores$Pseudotime, scores[,paste0("ME",i)], spar = 1), col = i)
          colors <- c(colors, i) 
        }
      }
  }
  legend("topleft", col = colors, legend = colors, lwd=.5, ncol = floor(length(colors)/3), cex = .75,lty=1)
}


## ------------------------------------------------------------------------
## ME boxplot
BoxMEByCluster <- function(srt, datExpr, net, which) {
  df <- data.frame(net$MEs[,which])
  rownames(df) <- rownames(datExpr)

  df$cluster <- srt@ident[rownames(df)]
  colnames(df) <- c("module", "cluster")

  ggplot(df, aes(cluster, module)) + geom_boxplot(aes(fill = cluster))
}

## ------------------------------------------------------------------------
## Uses a 2 mixture model to find cells that are + and - for some module
SelectMECells <- function(datExpr, net, which) {
  MEs <- net$MEs[, paste0("ME", which)]
  mm <- Mclust(MEs, G = 2)
  
  positive.clust <- which(mm$parameters$mean == max(mm$parameters$mean))
  negative.clust <- which(mm$parameters$mean == min(mm$parameters$mean))
  
  ans <- mm$z[,positive.clust]
  names(ans) <- rownames(datExpr)
  
  ans
}

## ------------------------------------------------------------------------
## Adds Module Eigengene Calculation to Seurat metadata
AddMEsToMetaData <- function(srt, net) {
  modules <- colnames(net$MEs)
  for(m in modules) {
    eigengenes <- net$MEs[, m]
    names(eigengenes) <- rownames(srt@meta.data)
    srt <- AddMetaData(srt, eigengenes, m)
  }
  srt
}

