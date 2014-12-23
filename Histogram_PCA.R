#' @title PCA on set of densities
#'
#' @description The function takes in a list of densities, runs a 
#' PCA on these densities, with common support.
#' @param dens list of density estimates
#' @param ... Arguments to be passed to \code{\link{prcomp}}.
#' @export
#' @return List of PCA, support, 
#' @examples 
#' N = 1000
#' n_i = 100
#' mix = 0.5
#' mean2= 1.5
#' sd2 = 2
#' mat1 = matrix(rnorm(N * mix * n_i), nrow=n_i)
#' mat2 = matrix(rnorm(N * (1-mix) * n_i, mean=mean2, sd = sd2), nrow=n_i)
#' mat = cbind(mat1, mat2)
#' dens = apply(mat, 2, density)
#' 
#' mat1 = matrix(runif(N * mix * n_i), nrow=n_i)
#' mat2 = matrix(runif(N * (1-mix) * n_i, max=mean2), nrow=n_i)
#' mat = cbind(mat1, mat2)
#' dens = apply(mat, 2, density) 
#' 
#' out = dens_pca(dens)
#' support = out$support
#' pr = out$pca
#' all.dens = out$all.dens
#' pct = out$pct_var
#' pcs = pr$x
#' ######################################
#' # Plotting over common support
#' ######################################
#' plot(support, rep(0, length=length(support)), 
#'      xlim = range(support), ylim= c(0, 1), type = "n")
#' apply(all.dens, 2, lines, x=support)
#' ######################################
#' # Plotting Variance explained for first 20 pcs
#' ######################################
#' plot(pct[1:20], type=c('l'))
#' points(pct[1:20], pch = 16, cex=0.5)
#' 
#' plot(x = support, y=pcs[,1], type='l')
#' lines(x = support, y=pcs[, 2], col="red")
#' 
dens_pca <- function(dens, # list of density estimates
  ... # Arguments to be passed to \code{\link{prcomp}}.
                        ){
    grids = sapply(dens, function(d){
      u = unique(round(diff(d$x), 8))
      stopifnot(length(u) == 1)
      u
    })
  ### get minimum grid for interpolation
  min.grid = min(grids)
  ####### get all ranges
  ranges = sapply(dens, function(d){
    range(d$x, na.rm=TRUE)
  })
  
  ### get total range of all supports
  all.range = c(min(ranges[1,]), max(ranges[2,]))
  
  ### need to get a little above the maximum
  real.max = ceiling(all.range[2]/min.grid) * min.grid
  real.max = real.max 
  
  ##### fulll support created
  support = seq(all.range[1], real.max, by=min.grid)
  
  ########## Approximating the density on the support
  all.dens = sapply(dens, function(d){
    y = approx(d$x, d$y, xout=support)$y
    y[is.na(y)] = 0
    y
  })
  
  #############################
  # All dens should be positive
  stopifnot(any(all.dens > 0))
  ### Run PCA
  pr = prcomp(all.dens, ...)
  
  pc1 = pr$x[,1]
  m.pc1 = max(abs(pc1))
  ind = which(pc1 == m.pc1 | pc1 == -m.pc1)
  #########################
  # Max of PC1 should be positive
  #########################
  if (length(ind) ==1 ){
    pcsign = sign(pc1[ind]) == -1
    if (pcsign){
      pr$x = -pr$x
      pr$rotation = -pr$rotation
    }
  }
  #### Get Percent variance explained
  pct = pr$sdev^2 / sum(pr$sdev^2)
#
# have.support = rowSums(all.dens) > 0
# sup.dens = all.dens[have.support, ]
# sup.sup = support[have.support]
# library(NMF)
# res <- nmf(sup.dens, 3)
#   
  return(list(pca = pr, 
              all.dens = all.dens, 
              support= support,
              pct_var = pct))
}


