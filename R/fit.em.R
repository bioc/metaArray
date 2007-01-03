"fit.em" <-
function(x, cl, threshold=1e-6){
  sup <- sum(cl==0) > 0 && sum(cl==1) > 0
  x <- as.numeric(x)
  len <- length(x)
  z <- rep(0,length(x))
  log.lik <- 1000
  lik.rec <- NULL
  num.iter <- 0
  a <- min(x,na.rm=TRUE); b <- max(x,na.rm=TRUE)
  err <- err.old <- 1
  if(!sup) {
    z <- find.init(x)
    Pi <- mean(z)
    mu <- sum((1-z)*x) / sum(1-z)
    sigmasq <- sum((1-z)*(x-mu)^2) / sum(1-z)
    tt <- len
    while(err > threshold) {
      num.iter <- num.iter + 1
      log.lik.old <- log.lik
      err.old <- err

      ## E Step
      est.u <- dunif(x,a,b)
      est.u.p <- Pi * est.u
      est.n <- dnorm(x,mu,sqrt(sigmasq))
      est.n.p <- (1-Pi) * est.n
      z <- est.u.p / (est.n.p + est.u.p)

      if(any(is.na(z))) stop("NA occurred in imputation\n")
      ## M Step
      mu <- sum((1-z)*x) / sum(1-z)
      sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)
      Pi <- sum(z) / len    
      sgn.z <- ifelse(x < mu, -1, 1)

      ## Likelihood
      est.u <- dunif(x,a,b)
      est.u.p <- Pi * est.u
      est.n <- dnorm(x,mu,sqrt(sigmasq))
      est.n.p <- (1-Pi) * est.n   
      log.lik <- sum(log(est.u.p + est.n.p))
      err <- abs(log.lik.old - log.lik)
      if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
    }    
  }
  else {
    tt <- sum(cl==0)
    z[cl==0] <- runif(tt,0,1)
    Pi <- mean(z[cl==0])
    mu <- sum((1-z)*x) / sum(1-z)
    sigmasq <- sum((1-z)*(x-mu)^2) / sum(1-z)

    while(err > threshold) {
      num.iter <- num.iter + 1
      log.lik.old <- log.lik
      err.old <- err

      ## E Step
      est.u <- dunif(x,a,b)
      est.u.p <- Pi * est.u
      est.n <- dnorm(x,mu,sqrt(sigmasq))
      est.n.p <- (1-Pi) * est.n
      z <- est.u.p / (est.n.p + est.u.p)
      z[cl==1] <- 0
      if(any(is.na(z))) stop("NA occurred in imputation\n")
      ## M Step
      mu <- sum((1-z)*x) / sum(1-z)
      sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)
      Pi <- mean(z[cl==0])    
      sgn.z <- ifelse(x < mu, -1, 1)

      ## Likelihood
      est.u <- dunif(x,a,b)
      est.u.p <- Pi * est.u
      est.n <- dnorm(x,mu,sqrt(sigmasq))
      est.n.p <- (1-Pi) * est.n   
      log.lik <- sum(log(est.u.p[cl==0] + est.n.p[cl==0])) + sum(log(est.n[cl==1]))
      err <- abs(log.lik.old - log.lik)
      if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
    }    
  }
  est.u.p <- Pi * dunif(x,a,b) 
  est.n.p <- (1-Pi) * dnorm(x,mu,sqrt(sigmasq))
  est.u.mu <- Pi * dunif(mu,a,b)
  est.n.mu <- (1-Pi) * dnorm(mu,mu,sqrt(sigmasq))
  z0 <- est.u.p / (est.n.p + est.u.p)
  zmu <- est.u.mu / (est.u.mu + est.n.mu)
  sgn.z0 <- ifelse(x < mu, -1, 1) 
  #loc <- (max(lik.rec) != lik.rec[length(lik.rec)])
  expr <- rep(0, len)
  expr <- (z0 - zmu) * sgn.z0
  return(list(expr=expr, a=a, b=b, sigmasq=sigmasq, mu=mu, Pi=Pi, lik.rec=lik.rec))
}

