## ------------------------------------------------------------------------
# Saves all feature plots
FeatureAll <- function(srt, marker.genes, out.dir) {
  for(i in marker.genes) {
    if(i %in% rownames(srt@scale.data)){
      print(FeaturePlot(srt, i, cols.use = c("darkblue", "red"), no.axes=T, no.legend = F))
      dev.print(pdf, file = sprintf(paste0(out.dir, "/%s.pdf"), i), width = 16, height = 9)
    }
  }
}

# Plots the top n DE genes (default n = 10)
FeatureTopDiff <- function(srt, mark, out.dir, n = 10) {
  mark %>% group_by(cluster) %>% top_n(n, avg_logFC) -> topn
  for(i in unique(topn$cluster)) {
      genes <- topn[topn$cluster == i,]$gene
      FeatureAll(srt, genes, sprintf(paste0(out.dir, "/cluster_%s"),i))
  }
}

# Saves all violin plots
ViolinAll <- function(srt, marker.genes, out.dir) {
  for(i in marker.genes) {
    if(i %in% rownames(srt@scale.data)){
      print(VlnPlot(srt, i))
      dev.print(pdf, file = sprintf(paste0(out.dir, "/%s.pdf"), i), width = 16, height = 5)
    }
  }
}

# Save all PCs
PCAll <- function(srt, PCs, out.path) {
  print(out.path)
  for(i in PCs) {
    print(PCHeatmap(srt, pc.use = i, cells.use = 100, do.balanced = T, label.columns = F))
    dev.print(pdf, file = sprintf(paste0(out.path, "/PC%s.pdf"), i), width = 9, height = 15)
  }
  
  # Saves jackstraw as well
  print(JackStrawPlot(srt, PCs = 1:ncol(srt@dr$pca@cell.embeddings)))
  dev.print(pdf, file = paste0(out.path, "/JackStraw.pdf"), width = 9, height = 25)
}

## ------------------------------------------------------------------------
# Heatmap of top DE genes based on average difference
# Default is n = 20 PCs with p < 0.01
diff.heatmap <- function(srt, mark, n = 20, p.cutoff = 0.01) {
  mark <- mark[mark$p_val_adj < p.cutoff, ]
  mark %>% group_by(cluster) %>% top_n(n, avg_logFC) -> topn
  
  
  
  # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
  DoHeatmap(srt, genes.use = topn$gene,  
            slim.col.label = TRUE, remove.key = TRUE, draw.line = T)
}

## ------------------------------------------------------------------------
## Plots a bunch of figures to a directory
FigureBomb <- function(srt, diffexp, markers, out.dir, do.pcs = T, do.tsne = T, 
                       do.markers = T, do.topclusts = T,
                       do.top5 = T, do.top20 = T) {
    # PCS
    if(do.pcs) {
        PCAll(srt, out.path = paste0(out.dir,"/Principal_Components"))
    }
    
    # TSNE
    if(do.tsne) {
        print(TSNEPlot(srt, do.label=T))
        dev.print(pdf, file = paste0(out.dir, "/TSNE.pdf"), width = 16, height = 9)
    }
    
    # Top 20 diff
    if(do.top20) {
        print(diff.heatmap(srt, diffexp, n = 20, p.cutoff = 1e-2))
        dev.print(pdf, file = paste0(out.dir, "/Top20_DiffExpr.pdf"), width = 16, height = 70)
    }
    
    # Top 5 diff
    if(do.top5) {
      print(diff.heatmap(srt, diffexp, n = 5, p.cutoff = 1e-2))
      
        dev.print(pdf, file = paste0(out.dir, "/Top5_DiffExpr.pdf"), width = 16, height = 30)
    }
    
    # Feature plots of marker genes
    if(do.markers){
        FeatureAll(srt, markers, paste0(out.dir, "/Feature_Plots/markers"))
    }
    
    # Features of top diff
    if(do.topclusts) {
        FeatureTopDiff(srt, diffexp, paste0(out.dir, "/Feature_Plots/clusters"), n = 10)
    }
    # Diff Expr table
    write.table(diffexp, file = paste0(out.dir, "/DiffExpr.tsv"), quote=F,sep="\t")
}

## ------------------------------------------------------------------------
## Plot Feature from metadata
TSNEPlotFeature <- function(srt, feats) {
    plotlist <- list()
    for(which in feats) {
        df <- data.frame(x = srt@dr$tsne@cell.embeddings[,1], 
                         y = srt@dr$tsne@cell.embeddings[,2], 
                         feat = as.numeric(srt@meta.data[, which]))
    
        plotlist[[which]] <- ggplot(df, aes(x = x, y = y, color = feat)) + 
            geom_point(size = 1) + 
            scale_colour_gradient(low = "lightgray", high = "red", na.value = "red") + 
            ggtitle(which)
    }
    
    multiplot(plotlist = plotlist, cols = floor(sqrt(length(plotlist))))
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## ------------------------------------------------------------------------
## Plots only markers that are unique to each cluster
PlotUnique <- function(srt, diffexp, out.dir, do.violin = F) {
    genes <- diffexp$gene
    unique.genes <- genes[!(duplicated(genes) | duplicated(genes, fromLast=T))]
    for(i in unique(diffexp$cluster)) {
        genes.plot <- unique(diffexp$gene[diffexp$cluster == i])
        genes.plot <- genes.plot[genes.plot %in% unique.genes]
        genes.plot <- genes.plot[1:(min(50, length(genes.plot)))]
        if(!do.violin) {
            FeatureAll(srt, genes.plot, paste0(out.dir, "/cluster_",i))
        } else {
            ViolinAll(srt, genes.plot, paste0(out.dir, "/cluster_",i))
        }
    }
}

## ------------------------------------------------------------------------
## Adds iteration to tSNE
RunMoreTSNE <- function(srt, n.iter, PCs, verbose = T) {
    res <- Rtsne(srt@dr$pca@cell.embeddings[, PCs], 
                 max_iter = n.iter, 
                 pca = F, 
                 Y_init = srt@dr$tsne@cell.embeddings,
                 verbose = verbose)$Y
    
    colnames(res) <- c("tSNE_1", "tSNE_2")
    rownames(res) <- rownames(srt@dr$pca@cell.embeddings)
    
    srt@dr$tsne@cell.embeddings <- res
    
    srt
}

## ------------------------------------------------------------------------
## Force directed layout plotting
PlotFDL <- function(srt, colors) {
  cors <- cor(as.matrix(srt@data[srt@var.genes,]))
  cors[cors == 1] <- 0
  thresh <- min(apply(cors, 2, max))
  cors <- cors >= thresh
  
  cl <- srt@ident
  names(cl) <- rownames(srt@meta.data)
  
  gg <- graph_from_adjacency_matrix(cors, mode="undirected")
  gg <- simplify(gg)
  gg <- delete.vertices(gg, which(degree(gg) == 0))
  V(gg)$color <- "#2a5e8c"
  
  print("Plotting...")
  plot(gg, vertex.size = 2, vertex.label=NA)
  legend('topleft',legend=paste0("Cluster ", 1:length(unique(cl))), col=1:length(unique(cl)), pt.bg = 1:length(unique(cl)), pt.cex = 1, pch = 19)
}

## ------------------------------------------------------------------------
## Calculates distance between centroids of clusters in PC space
centroiddist <- function(srt, PCs){
  idents <- unique(srt@ident)
  ans <- NULL
  for(i in idents) {
    ans <- rbind(ans, colMeans(srt@dr$pca@cell.embeddings[srt@ident == i,PCs]))
  }
  rownames(ans) <- idents
  as.matrix(dist(ans))
}


