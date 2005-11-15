require(MergeMaid)

intcor <- function(merged) {     
  l <- length(merged)
  if(l<=1) stop("More than a single study is required to run integrative correlation\n")
  inter <- intersection(merged)
  common.gname <- rownames(exprs(inter))
  expr <- cl <- NULL
  for(i in 1:l) {
    expr[[i]] <- exprs(merged[i])[match(common.gname,rownames(exprs(merged[i]))),]
    cl[[i]] <- rep(0,dim(expr[[i]])[2])
  }
  nr <- length(common.gname)
  nc <- rep(0,l)
  for(i in 1:l) nc[i] <- length(cl[[i]])
  
  res <-
    .C("intcor",as.double(unlist(expr)),
       as.integer(unlist(cl)),
       as.integer(l),
       as.integer(nr),
       as.integer(nc),
       correl=double(nr),
       paircor=double((l*(l-1)/2)*nr),
       PACKAGE="metaArray"
      )
  tmp.paircor <- matrix(res$paircor, nrow=nr)
  rownames(tmp.paircor) <- common.gname
  obj <- list(pair.cor=tmp.paircor, avg.cor=res$correl)
  return(obj)
}



  
