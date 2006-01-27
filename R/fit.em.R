"fit.em" <-
function(x, cl, threshold=0.00001){ 
  conv <- T
  dec <- F
  abnormal <- F
  x <- as.numeric(x)
  nx <- length(x)
  ind1 <- (1:nx)[!is.na(x) & cl==1]
  ind2 <- (1:nx)[!is.na(x) & cl==0]
  x2 <- x[ind1]
  tt <- length(x2)
  z <- rep(0,length(x2))
  id <- find.init(x2)
  log.lik <- 1000
  lik.rec <- NULL
  num.iter <- 0
  if(all(id==1) | all(id==0)) {
    conv <- F
    abnormal <- T
    lik.rec <- 0
  }
  z[id == 1] <- 1
  a <- min(x,na.rm=T); b <- max(x,na.rm=T)
  a <- a - (b-a)*.0001; b <- b + (b-a)*.0001  
  err <- err.old <- 1
  pi <- sum(z) / tt 
  mu <- sum((1-z)*x2) / sum(1-z)
  sigmasq <- sum((1-z)*(x2-mu)^2) / sum(1-z)
  while(err > threshold & conv == T) {
    num.iter <- num.iter + 1
    log.lik.old <- log.lik
    err.old <- err
    z.old <- z
   ## E Step
    est.u <- dunif(x2,a,b)
    est.u.p <- pi * est.u
    est.n <- dnorm(x2,mu,sqrt(sigmasq))
    est.n.p <- (1-pi) * est.n
    z <- est.u.p / (est.n.p + est.u.p)
    if(any(is.na(z))) conv <- F
   ## M Step
    pi <- sum(z) / tt 
    sigmasq <- sum((1-z)*((x2-mu)^2)) / sum(1-z)
    mu <- sum((1-z)*x2) / sum(1-z)
    sgn.z <- ifelse(x2 < mu, -1, 1)  
    est.u[est.u==0] <- 10e-2 
    est.n[est.n==0] <- 10e-2
    log.lik <- sum(log(est.u.p + est.n.p))
    err <- abs(log.lik.old - log.lik)
    if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
    dec <- ifelse(log.lik.old > log.lik, T, F)
    if(any(pi==0 | pi==1)) err <- 0
  }
  x3 <- x[ind2]
  tt <- length(x3)
  z0 <- rep(0,tt)
  est.u.p <- pi * dunif(x3,a,b)
  est.n.p <- (1-pi) * dnorm(x3,mu,sqrt(sigmasq))
  z0 <- est.u.p / (est.n.p + est.u.p)
  sgn.z0 <- ifelse(x3 < mu, -1, 1) 
  loc <- (max(lik.rec) != lik.rec[length(lik.rec)])
  if(!conv | abnormal) {
    expr <- rep(NA,nx)
    mu <- sigmasq <- pi <- NA
    loc <- F
    lik.rec <- 0
  }
  else {
    expr <- rep(NA, nx)
    expr[ind1] <- (z - min(z)) * sgn.z
    expr[ind2] <- (z0 - min(z)) * sgn.z0  
  }
  return(list(expr=expr, a=a, b=b, sigmasq=sigmasq, mu=mu, pi=pi, conv=conv, 
dec=dec, loc=loc, lik.rec=lik.rec, abnormal=abnormal))
}

