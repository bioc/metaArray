"find.init" <-
function(z, width = 1) {
  len <- length(z)
  init <- rep(0,len)
  for(i in 1:len) {
    test.z <- z[i]
    train.z <- z[-i]
    med.z <- median(train.z)
    sigmasq <- sum((train.z-med.z)^2)/(len-1)
    dr <- abs(test.z - med.z)
    dthe <- width * sqrt(sigmasq)
    init[i] <- as.numeric(dr > dthe)  
  }
  return(init)
}

