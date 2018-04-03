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
PCAll <- function(srt, out.path) {
  print(out.path)
  for(i in 1:40) {
    print(PCHeatmap(srt, pc.use = i, cells.use = 100, do.balanced = T, label.columns = F))
    dev.print(pdf, file = sprintf(paste0(out.path, "/PC%s.pdf"), i), width = 9, height = 15)
  }
  
  # Saves jackstraw as well
  print(JackStrawPlot(srt, PCs = 1:50))
  dev.print(pdf, file = paste0(out.path, "/JackStraw.pdf"), width = 9, height = 25)
}

## ------------------------------------------------------------------------
diff.heatmap <- function(srt, mark, n = 20, p.cutoff = 1e-5) {
  mark <- mark[mark$p_val_adj < p.cutoff, ]
  mark %>% group_by(cluster) %>% top_n(n, avg_logFC) -> topn
  
  
  
  # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
  DoHeatmap(srt, genes.use = topn$gene,  
            slim.col.label = TRUE, remove.key = TRUE, draw.line = T)
}

## ------------------------------------------------------------------------
FindClustersSemi <- function(srt, mtx, k.iter = 5:50, pc.use = 1:20) {
    markers <- rownames(mtx)
    nmark <- nrow(mtx)
    print(markers)
    ans <- NULL
    minscore <- 100000000
    which.min <- NULL
    for(k in k.iter) {
        # Find clusters
        print(paste0("k = ",k))
        tmp <- FindClusters(srt, reuse.SNN = F, dims.use = pc.use, k.param = k, print.output = F)
        nclust <- length(unique(tmp@ident))
        print(paste0("Num clusters = ", nclust))
        
        # Finds DE genes
        mark <- FindAllMarkers(tmp, only.pos=T, thresh.use = 1)
        
        # Gets top 10 DE genes for each cluster
        topdiff <- data.frame(mark %>% group_by(cluster) %>% top_n(100,avg_logFC))

        # Builds a 1-0 matrix saying if gene is marker for each cluster
        mark.clust <- matrix(NA, nrow = nmark, ncol = nclust)
        for(i in 1:nmark) {
            for(j in 1:nclust) {
                mark.clust[i,j] <- sum(markers[i] %in% topdiff[topdiff$cluster == j - 1,]$gene)
            }
        }
        
        rownames(mark.clust) <- markers
        colnames(mark.clust) <- 1:nclust
        print(mark.clust)
        # Finds how many times each pair of genes is in the matrix
        mark.dist <- matrix(NA, nrow = nmark, ncol = nmark)
        rownames(mark.dist) <- colnames(mark.dist) <- markers
        for(i in 1:nmark) {
            for(j in 1:nmark) {
                mark.dist[i,j] <- sum(mark.clust[i,] & mark.clust[j,])
            }
        }
        
        print(mtx)
        score <- sum((mtx[markers, markers] - mark.dist[markers, markers])^2)
        print(paste0("score = ",score))
        if(score < minscore) {
            
            print("--------NEW MIN!!!!--------")
            minscore <- score
            which.min <- k
            ans <- tmp
        }
    }
    
    print(paste0("min score = ", minscore))
    print(paste0("which k min = ", which.min))
    ans
}

## ------------------------------------------------------------------------
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
RunMoreTSNE <- function(srt, n.iter = 1000, pcs.use = 1:50, verbose = T) {
    res <- Rtsne(srt@dr$pca@cell.embeddings[, pcs.use], 
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
ViolinMarkers <- function(srt, markers, ncol = 1) {
  df <- NULL
  for(i in markers){
    vals <- data.frame(expr = srt@data[i,])
    vals$gene <- as.character(i)
    vals$cluster <- srt@ident
    
    df <- rbind(df, vals)
  }
  
  df$cluster <- factor(df$cluster)
  print(head(df))
  
  ggplot(df, aes(x = cluster, y = expr)) +
    geom_violin(aes(fill = cluster)) + facet_wrap(~gene, ncol = ncol, scales = "free")
}

## ------------------------------------------------------------------------
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
AddMito <- function(srt) {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = srt@data), value = TRUE)
  percent.mito <- Matrix::colSums(srt@raw.data[mito.genes, ])/Matrix::colSums(srt@raw.data)

  srt <- AddMetaData(object = srt, metadata = percent.mito, col.name = "percent.mito")
}

## ------------------------------------------------------------------------
ProjectSeurat <- function(sa, sb, pref.a = "Hu", pref.b = "Org") {
  # Subsets only genes present in both
  write("Subsetting genes...", stderr())
  genes.both <- intersect(rownames(sa@data), rownames(sb@data))
  sa@data <- sa@data[genes.both,]
  sb@data <- sb@data[genes.both,]
  sa@raw.data <- sa@data[genes.both,]
  sb@raw.data <- sb@data[genes.both,]
  
  # Prepends prefix to data, raw.data, meta.data and cell.embeddings
  write("Prepending names to stuff...", stderr())
  colnames(sa@data) <- paste0(pref.a, "_", colnames(sa@data))
  colnames(sb@data) <- paste0(pref.b, "_", colnames(sb@data))
  colnames(sa@raw.data) <- paste0(pref.a, "_", colnames(sa@raw.data))
  colnames(sb@raw.data) <- paste0(pref.b, "_", colnames(sb@raw.data))
  rownames(sa@meta.data) <- paste0(pref.a, "_", rownames(sa@meta.data))
  rownames(sb@meta.data) <- paste0(pref.b, "_", rownames(sb@meta.data))
  rownames(sa@dr$pca@cell.embeddings) <- paste0(pref.a, "_", rownames(sa@dr$pca@cell.embeddings))
  rownames(sb@dr$pca@cell.embeddings) <- paste0(pref.b, "_", rownames(sb@dr$pca@cell.embeddings))
  
  # Creates merged origin list 
  merged.origin <- c(rep(pref.a, nrow(sa@meta.data)), rep(pref.b, nrow(sb@meta.data)))
  names(merged.origin) <- c(rownames(sa@meta.data), rownames(sb@meta.data))
  
  # PCA projection
  write("Projecting Seurat B onto Seurat A", stderr())
  genes.pc <- intersect(genes.both, rownames(sa@dr$pca@gene.loadings))
  proj.mat <- sa@dr$pca@gene.loadings[genes.pc,]
  
  # Projection of Seurat B onto Seurat A PC Space
  sb.proj <- as.matrix(t(proj.mat) %*% sb@data[genes.pc, ])
  
  write("Adding all Seurat  B data onto Seurat A", stderr())
  sa@dr$pca@cell.embeddings <- rbind(sa@dr$pca@cell.embeddings, t(sb.proj))
  sa@data <- cbind(sa@data, sb@data)
  sa@raw.data <- cbind(sa@raw.data, sb@raw.data)

  write("Adding relevant metadata of Seurat B onto Seurat A", stderr())
  meta.both <- intersect(colnames(sa@meta.data), colnames(sb@meta.data))
  sa@meta.data <- rbind(sa@meta.data[, meta.both], sb@meta.data[, meta.both])
  
  write("Adding Merged Origin Metadata", stderr())
  sa <- AddMetaData(sa, merged.origin, "merged.origin")
  
  write("Overriding idents...", stderr())
  sa@ident <- factor(c(paste0(pref.a,"_", sa@ident), paste0(pref.b, "_", sb@ident)))
  names(sa@idnet) <- rownames(sa@meta.data)
  
  sa

}

## ------------------------------------------------------------------------
avgdist <- function(srt, PCs = 1:19){
  idents <- unique(srt@ident)
  ans <- matrix(NA, nrow = length(idents), ncol = length(idents))
  for(i in idents) {
    for(j in idents) {
      mmat <- matrix(NA, nrow=  sum(srt@ident ==i), ncol = sum(srt@ident == j))
      va <- srt@dr$pca@cell.embeddings[srt@ident == i, PCs]
      vb <- srt@dr$pca@cell.embeddings[srt@ident == j, PCs]
      for(k in 1:nrow(va)) {
        for(l in 1:nrow(vb)) {
          mmat[k,l] <- sum((va[k,] - vb[l,])^2)
        }
      }
      print(mmat)
      print(dim(mmat))
      ans[as.integer(i),as.integer(j)] <- mean(mmat)
    }
  }
  rownames(ans) <- colnames(ans) <- paste0("H17w_", idents)
  ans
}

centroiddist <- function(srt, PCs = 1:19){
  idents <- unique(srt@ident)
  ans <- NULL
  for(i in idents) {
    ans <- rbind(ans, colMeans(srt@dr$pca@cell.embeddings[srt@ident == i,PCs]))
  }
  rownames(ans) <- idents
  as.matrix(dist(ans))
}


## ------------------------------------------------------------------------


