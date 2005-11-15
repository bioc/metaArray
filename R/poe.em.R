"poe.em" <-
function(mat, cl=NULL, threshold=0.00001, every=100) {
  mat <- as.matrix(mat)
  nc <- ncol(mat)
  nr <- nrow(mat)
  if(is.null(cl)) cl <- rep(1,dim(mat)[2])
  cat("Number of Samples:", nc, "\n")
  cat("Number of Genes:", nr, "\n")
  cat("This model assumes that the samples are centered and scaled.\n")
  new.mat <- matrix(0,nr,nc)
  conv <- dec <- loc <- abnormal <- rep(F,nr)
  med.expr <- apply(mat,1,median,na.rm=T)
  new.mat <- sweep(mat,1,med.expr)
  for(i in 1:nr) {
    if(sum(is.na(as.numeric(mat[i,]))) > .25 * nc ) stop("More than 25% missing values for gene", i, "\n")
    zvec <- fit.em(as.numeric(mat[i,]), cl, threshold=threshold) 
    new.mat[i,] <- zvec$expr
    conv[i] <- zvec$conv  
    dec[i] <- zvec$dec
    loc[i] <- zvec$loc
    abnormal[i] <- zvec$abnormal
    if(i%%every==0) cat(i, "genes fitted\n")   
  }
  checkid <- (1:nr)[abnormal]
  len.id <- length(checkid)
  rownames(new.mat) <- rownames(mat)
  colnames(new.mat) <- colnames(mat)
  cat("The algorithm has shown evidences of\n")
  cat("local maxima for", sum(loc), "genes\n")
  cat("decreasing log-likelihood", sum(dec), "genes\n")
  if(len.id!=0) {
    cat("EM failed to find starting points for:\n")
    for(i in 1:len.id) cat(checkid[i], rownames(mat)[i], "\n")
  }
  return(list(data=new.mat, conv=conv, dec=dec, loc=loc))
}

