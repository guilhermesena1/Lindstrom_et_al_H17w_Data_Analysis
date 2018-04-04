## ------------------------------------------------------------------------
## Bhattacharyya distance
bd <- function(res, cl1, cl2) {
  sigma1 <- res$parameters$variance$sigma[,,cl1]
  sigma2 <- res$parameters$variance$sigma[,,cl2]
  sigma <- (sigma1 + sigma2)/2
  mu1 <- res$parameters$mean[,cl1]
  mu2 <- res$parameters$mean[,cl2]
  
  d <- as.numeric(det(sigma))
  d1 <- as.numeric(det(sigma1))
  d2 <- as.numeric(det(sigma2))
  
  t(mu1 - mu2) %*% solve(sigma)%*%(mu1 - mu2)/8 #+ log(d/sqrt(d1*d2))/2
}

## ------------------------------------------------------------------------
## KL Divergence between two GMM clusters
kl <- function(res, cl1, cl2) {
  sigma1 <- res$parameters$variance$sigma[,,cl1]
  sigma2 <- res$parameters$variance$sigma[,,cl2]
  mu1 <- res$parameters$mean[,cl1]
  mu2 <- res$parameters$mean[,cl2]
  
  (normdiff(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, method="KL") +
  normdiff(mu1 = mu2, sigma1 = sigma2, mu2 = mu1, sigma2 = sigma1, method="KL"))/2
}

## ------------------------------------------------------------------------
## Hellinger Distance Between two GMM Clusters
hellinger <- function(res, cl1, cl2) {
  sigma1 <- res$parameters$variance$sigma[,,cl1]
  sigma2 <- res$parameters$variance$sigma[,,cl2]
  mu1 <- res$parameters$mean[,cl1]
  mu2 <- res$parameters$mean[,cl2]
  
  normdiff(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, method="Hellinger")
}


## ------------------------------------------------------------------------
## BD of two Mclust MLE objects
bd.mle <- function(res1, res2) {
  sigma1 <- res1$parameters$variance$Sigma
  sigma2 <- res2$parameters$variance$Sigma
  sigma <- (sigma1 + sigma2)/2
  mu1 <- res1$parameters$mean
  mu2 <- res2$parameters$mean
  
  d <- as.numeric(det(sigma))
  d1 <- as.numeric(det(sigma1))
  d2 <- as.numeric(det(sigma2))
  
  t(mu1 - mu2) %*% solve(sigma)%*%(mu1 - mu2)/8 + log(d/sqrt(d1*d2))/2
}

## ------------------------------------------------------------------------
## Calculates Bhattacharyya Distance of two Seurat clusters by estimating Gaussian MLE
bd.seurat <- function(srt, PCs = 1:30) {
  idents <- unique(srt@ident)
  mcl <- list()
  for(id in idents) {
    mcl[[id]] <- Mclust(srt@dr$pca@cell.embeddings[srt@ident == id, PCs], G = 1)
  }
  ans <- matrix(0, nrow = length(idents), ncol = length(idents))
  rownames(ans) <- idents
  colnames(ans) <- idents
  
  for(i in 1:(length(idents) - 1)) {
    for(j in (i+1):length(idents)) {
      ans[idents[i], idents[j]] <- ans[idents[j], idents[i]] <- bd.mle(mcl[[idents[i]]], mcl[[idents[j]]])
    }
  }
  ans
}

## ------------------------------------------------------------------------
## Functions for pairwise distances
pairwise.bd <- function(res, prefix = "Whole_"){
  cls <- unique(res$classification)
  ans <- matrix(NA, nrow = length(cls), ncol = length(cls))
  rownames(ans) <- colnames(ans) <- 1:length(cls)
  
  for(i in cls) {
    for(j in cls) {
      ans[i,j] <- bd(res, i,j)
    }
  }
  
  ans
}

pairwise.kl <- function(res, prefix = "Whole_"){
  cls <- unique(res$classification)
  ans <- matrix(NA, nrow = length(cls), ncol = length(cls))
  rownames(ans) <- colnames(ans) <- 1:length(cls)
  
  for(i in cls) {
    for(j in cls) {
      ans[i,j] <- kl(res, i,j)
    }
  }
  
  ans
}

pairwise.hellinger <- function(res, prefix = "Whole_"){
  cls <- unique(res$classification)
  ans <- matrix(NA, nrow = length(cls), ncol = length(cls))
  rownames(ans) <- colnames(ans) <- 1:length(cls)
  
  for(i in cls) {
    for(j in cls) {
      ans[i,j] <- hellinger(res, i,j)
    }
  }
  
  ans
}


## ------------------------------------------------------------------------
## Plots hclust for pairwise bd
plot.pairwise.bd <- function(bds) {
  plot(hclust(as.dist(bds)), hang = -1, xlab = "Cluster", ylab = "Bhattacharyya distance", main = "GMM Cluster Hierarchy")
  
}

## Plots MST with BD
plot.bd.mst <- function(bds) {
  gg <- simplify(graph_from_adjacency_matrix(bds, mode="undirected", weighted = T))
  gg <- mst(gg)
  plot(gg, vertex.size = 4, vertex.label.cex = .8, vertex.label.dist = .3, 
       vertex.color="red", vertex.frame.color="black")
}


