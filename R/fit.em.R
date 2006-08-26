"fit.em" <-
function(x, cl=NULL, threshold=1e-10){ 
  conv <- TRUE
  dec <- FALSE
  abnormal <- FALSE
  x <- as.numeric(x)
  tt <- length(x)
  z <- rep(0,length(x))
  if(is.null(cl)) cl <- rep(0, length(x))
  if(all(cl==0)) id <- find.init(x)
  else id <- cl    
  log.lik <- 1000
  lik.rec <- NULL
  num.iter <- 0
  if(all(id==1) | all(id==0)) {
    conv <- FALSE
    abnormal <- TRUE
    lik.rec <- 0
  }
  if(all(cl==0)) z[id == 1] <- 1
  else z <- rep(0.5,tt)

  a <- min(x,na.rm=TRUE); b <- max(x,na.rm=TRUE)
  #a <- a - (b-a)*.001; b <- b + (b-a)*.001  
  err <- err.old <- 1

  if(all(cl==0)) pi <- sum(z) / tt
  else pi <- sum(cl) / tt

  mu <- sum((1-z)*x) / sum(1-z)
  sigmasq <- sum((1-z)*(x-mu)^2) / sum(1-z)

  while(err > threshold & conv == TRUE) {
    num.iter <- num.iter + 1
    log.lik.old <- log.lik
    err.old <- err
    z.old <- z
   ## E Step
    est.u <- dunif(x,a,b)
    est.u.p <- pi * est.u
    est.n <- dnorm(x,mu,sqrt(sigmasq))
    est.n.p <- (1-pi) * est.n
    z <- est.u.p / (est.n.p + est.u.p)

    if(any(is.na(z))) conv <- FALSE
   ## M Step
    if(any(cl==1)) z[cl==0] <- 0
    pi <- sum(z) / tt
    if(any(cl==1)) pi <- max(pi, 1 / tt)
    #if(any(cl==1)) pi <- sum(z) / sum(cl)
    mu <- sum((1-z)*x) / sum(1-z)
    sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)

    sgn.z <- ifelse(x < mu, -1, 1)  
    #est.u[est.u==0] <- 1e-15 
    #est.n[est.n==0] <- 1e-15

    est.u <- dunif(x,a,b)
    est.u.p <- pi * est.u
    est.n <- dnorm(x,mu,sqrt(sigmasq))
    est.n.p <- (1-pi) * est.n
    
    log.lik <- sum(log(est.u.p + est.n.p))
    err <- abs(log.lik.old - log.lik)
    if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
    dec <- ifelse(log.lik.old > log.lik, TRUE, FALSE)
    if(any(pi==0 | pi==1)) err <- 0
  }

  if(any(cl==1)) pi <- max(1/tt,pi)
  est.u.p <- pi * dunif(x,a,b) 
  est.n.p <- (1-pi) * dnorm(x,mu,sqrt(sigmasq))
  est.u.mu <- pi * dunif(mu,a,b)
  est.n.mu <- (1-pi) * dnorm(mu,mu,sqrt(sigmasq))
  z0 <- est.u.p / (est.n.p + est.u.p)
  zmu <- est.u.mu / (est.n.mu + est.u.mu)
  sgn.z0 <- ifelse(x < mu, -1, 1) 
  loc <- (max(lik.rec) != lik.rec[length(lik.rec)])
  if(!conv | abnormal) {
    expr <- rep(NA,tt)
    mu <- sigmasq <- pi <- NA
    loc <- FALSE
    lik.rec <- 0
  }
  else {
    expr <- rep(NA, tt)
    expr <- (z0-zmu) * sgn.z0
  }
  return(list(expr=expr, a=a, b=b, sigmasq=sigmasq, mu=mu, pi=pi, conv=conv, dec=dec, loc=loc, lik.rec=lik.rec, abnormal=abnormal))
}

