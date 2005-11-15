"Zscore" <-
function(merged, pheno=NULL, permute=0, verbose=TRUE) {     
  l <- length(merged)
  inter <- intersection(merged)
  common.gname <- rownames(exprs(inter))
  expr <- cl <- NULL
  if(is.null(pheno)) {
    cat("Pheno data is assumed to be in the first column of phenoData slot\n")
    pheno <- rep(1,l)
  }
  for(i in 1:l) {
    expr[[i]] <- exprs(merged[i])[match(common.gname,rownames(exprs(merged[i]))),]
    cl[[i]] <- pData(merged[i])[,pheno[i]]
  }

  # Logistics: number of samples, class label
  n <- length(expr)
  for(i in 1:n) {
    if(length(unique(cl[[i]]))!=2) {
      cat("Error in labels from data", i, "\n") 
      stop("Class label allows only two unique elements\n")
    }
  }
  mu0 <- mu1 <- sigma0 <- sigma1 <- NULL
  lab <- unique(cl[[1]]); lab <- lab[order(lab)]
  cat(lab[1], "marked as 0\n")
  cat(lab[2], "marked as 1\n")
  cat("Contrast will be", lab[2], "-", lab[1], "\n")
  for(i in 2:n) {
    tmp <- unique(cl[[i]]); tmp <- tmp[order(tmp)]
    if(sum(lab == tmp)!=2)
      stop("Disparate class labels between datasets\n")
  }

  # Change class labels to 0-1
  for(i in 1:n) {
    cl[[i]][cl[[i]]==lab[1]] <- 0
    cl[[i]][cl[[i]]==lab[2]] <- 1
  }
  
  # calculate mu and sigma's within each class label
  # for all datasets. make a list object for it.
  # Script smoothed mean-variance curve and penalized var estimates

  # Calculate the overall contrast by formula.
  nr <- length(common.gname)
  nc <- rep(0,l)
  for(i in 1:l) nc[i] <- length(cl[[i]])
  
  res <-
    .C("contr",as.double(unlist(expr)),
       as.integer(unlist(cl)),
       as.integer(l),
       as.integer(nr),
       as.integer(nc),
       as.integer(permute),
       z=double(nr),
       p=double(nr),
       PACKAGE="metaArray"
      )
  obj <- as.data.frame(cbind(Zscore=res$z,Pvalue=res$p))
  rownames(obj) <- common.gname
  return(obj)
}

