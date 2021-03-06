## ------------------------------------------------------------------------
## Converts Seurat Object into Monocle Object and Run Monocle Pipeline
monocle.from.seurat <- function(srt){
  cerr("Building Monocle Object...")
  
  # Annotation
  srt@meta.data$ident <- srt@ident
  anno <- new("AnnotatedDataFrame", data = srt@meta.data)
  mon <- newCellDataSet(as.matrix(srt@raw.data), phenoData = anno)
  
  
  cerr("Normalizing...")
  mon <- setOrderingFilter(mon, srt@var.genes)
  mon <- estimateSizeFactors(mon)
  mon <- estimateDispersions(mon)

  cerr("Running ICA...")
  mon <- reduceDimension(mon, max_components = 2)
  
  mon <- orderCells(mon)
  mon
}

## ------------------------------------------------------------------------
## Filters differential gene test by n cells cutoff and q val cutoff
do.dgt <- function(srt, mon) {
  genes.use <- srt@var.genes[srt@var.genes %in% rownames(mon@assayData$exprs)]
  differentialGeneTest(mon[genes.use,], relative_expr = F)
}

dgt.filter <- function(dgt, qval.cutoff = 1, ncells.cutoff = 0) {
  dgt <- dgt[dgt$qval < qval.cutoff, ]
  dgt <- dgt[dgt$num_cells_expressed > ncells.cutoff, ]
  dgt <- dgt[order(dgt$qval),]
  
  dgt
}

## ------------------------------------------------------------------------
## Subpath analysis: Selects only a subset of states from a monocle object
mon.sub <- function(mon, states) {
    cerr("Subsetting...")
    cells.keep <- rownames(mon@phenoData@data)[mon@phenoData@data$State %in% states]

    ret <- mon[, cells.keep]
    cerr("Doing ICA...")
    ret <- reduceDimension(ret, max_components = 2)
    cerr("Ordering...")
    ret <- orderCells(ret)
    ret
}

## ------------------------------------------------------------------------
## Plots expression heatmap for top genes using pheatmap
monocle.heatmap <- function(srt, mon, diff) {
  # gets cells in order
  cells.by.pseudotime <- colnames(srt@data)[colnames(srt@data) %in% colnames(mon@assayData$exprs)]
  cells.by.pseudotime <- cells.by.pseudotime[order(mon@phenoData@data[cells.by.pseudotime,"Pseudotime"])]
  genes.top <- rownames(diff)
  anno <- data.frame(Pseudotime = mon@phenoData@data[cells.by.pseudotime, "Pseudotime"])
  rownames(anno) <- cells.by.pseudotime
  
  #to.plot <- as.matrix(srt@data[genes.top, cells.by.pseudotime])
  #for(i in 1:nrow(to.plot)) {
  #  to.plot[i,] <- scale(to.plot[i,])
  #}
  to.plot[to.plot > 5] <- 5
  to.plot[to.plot < -5] <- -5
  
  
  pheatmap(to.plot, 
           cluster_cols = F, cluster_rows = T, 
           show_colnames = F, show_rownames = T, 
           annotation = anno)
}

## ------------------------------------------------------------------------
## Plots and saves pseudotime heatmaps
PseudotimeFeatureAll <- function(mon, srt, genes, out.dir) {
  for(i in genes) {
    print(PseudotimeFeaturePlot(mon, srt, i))
    dev.print(pdf, file = sprintf(paste0(out.dir,"/%s.pdf"), i), width = 16, height = 9)
  }
}

## ------------------------------------------------------------------------
GeneAll <- function(mon, genes, out.dir) {
  for(i in genes) {
    print(plot_genes_in_pseudotime(mon[i,], relative_expr = F, color_by="tree.ident"))
    dev.print(pdf, file = sprintf(paste0(out.dir, "/%s.pdf"), i))
  }
}


## ------------------------------------------------------------------------
PseudotimeFeaturePlot <- function(mon, srt, gene) {
  ica <-t(reducedDimS(mon))
  df <- data.frame(ica1 = ica[,1], ica2 = ica[,2], expr = srt@data[gene, rownames(ica)])
  
  ggplot(aes(x = ica1, y = ica2, colour = expr), data = df) + 
    geom_point() + 
    scale_color_gradient(low="darkblue", high="red") +
    ggtitle(gene)
}

## ------------------------------------------------------------------------
## Plots monocle tree when reduced to 3D
Plot3DCellTrajectory <- function(mon, srt, color_by = "ident", color.pal = NULL) {
  dims <- t(reducedDimS(mon))
  col <- NULL
  cells <- rownames(pData(mon))
  # Color by state
  if(color_by == "State") {
    col <- factor(pData(mon)$State)
    
  # Color by Seurat cluster
  } else if(color_by == "ident") {
    col <- factor(srt@ident[cells])
  
  # Color by gene expression
  } else if(color_by %in% rownames(srt@data)){
    col <- srt@data[color_by,cells]
  } else {
    stop(paste0("Dont know how to parse color_by: ", color_by))
  }
  
  # Makes data frame
  df <- data.frame(DDR_1 = dims[,1], DDR_2 = dims[,2], DDR_3 = dims[,3], col = col)
  
  # If no color palette, use default plotly
  if(is.null(color.pal)) {
    p <- plot_ly(data = df, x = ~DDR_1, y = ~DDR_2, z = ~DDR_3, 
          text = col, 
          color = ~col,
          marker = list(size = 2)
          )
    
  # Else use color palette as inscrution
  } else {
     p <- plot_ly(data = df, x = ~DDR_1, y = ~DDR_2, z = ~DDR_3, 
          text = col, 
          marker = list(size = 2, color = color.pal[as.integer(col)])
          )
  }
  p <- p %>% layout(title = color_by)
  p
}

## ------------------------------------------------------------------------
## Plots a list of genes in 3D
Plot3DAll <- function(mon, srt, genes, out.dir) {
  for(gene in genes) {
    p <- Plot3DCellTrajectory(mon, srt, gene)
    htmlwidgets::saveWidget(as_widget(p), paste0(out.dir, "/", gene, ".html"))
  }
}

## ------------------------------------------------------------------------
## Plots top DE genes in 3D Monocle
Pseudotime3DTopDiff <- function(mon, srt, mark, out.dir, n) {
  mark %>% group_by(cluster) %>% top_n(n, avg_logFC) -> topn
  for(i in unique(topn$cluster)) {
      genes <- topn[topn$cluster == i,]$gene
      Plot3DAll(mon, srt, genes, paste0(out.dir, "/cluster_", i))
  }
}

