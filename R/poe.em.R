"poe.em" <-
function(mat, cl, threshold=0.00001, every=100) {
  mat <- as.matrix(mat)
  nc <- ncol(mat)
  nr <- nrow(mat)
  if(all(is.null(cl))) cl <- rep(0,dim(mat)[2])
  cat("Number of Samples:", nc, "\n")
  cat("Number of Genes:", nr, "\n")
  cat("This model assumes that the samples are centered and scaled.\n")
  new.mat <- matrix(0,nr,nc)
  med.expr <- apply(mat,1,median,na.rm=TRUE)
  new.mat <- sweep(mat,1,med.expr)
  for(i in 1:nr) {
    if(sum(is.na(as.numeric(mat[i,]))) > .25 * nc ) stop("More than 25% missing values for gene", i, "\n")
    zvec <- fit.em(as.numeric(mat[i,]), cl, threshold=threshold) 
    new.mat[i,] <- zvec$expr
    if(i%%every==0) cat(i, "genes fitted\n")   
  }
  rownames(new.mat) <- rownames(mat)
  colnames(new.mat) <- colnames(mat)
  return(list(data=new.mat))
}

