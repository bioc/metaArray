"fit.em" <-
function(x, cl=NULL, threshold=0.0001){ 
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
  else z <- cl/2

  a <- min(x,na.rm=TRUE); b <- max(x,na.rm=TRUE)
  #a <- a - (b-a)*.001; b <- b + (b-a)*.001  
  err <- err.old <- 1

  if(all(cl==0)) pi <- sum(z) / tt
  else pi <- sum(z[cl==1]) / sum(cl)

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
    if(all(cl==0)) z <- est.u.p / (est.n.p + est.u.p)
    else {
      z <- est.u.p / (est.n.p + est.u.p) #[cl==1]
      z[cl==0] <- 0
    }

    if(any(is.na(z))) conv <- FALSE
   ## M Step
    if(all(cl==0)) pi <- sum(z) / tt
    else pi <- sum(z[cl==1]) / sum(cl)
    mu <- sum((1-z)*x) / sum(1-z)
    sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)

    sgn.z <- ifelse(x < mu, -1, 1)  
    est.u[est.u==0] <- 10e-2 
    est.n[est.n==0] <- 10e-2
    
    if(all(cl==0)) log.lik <- sum(log(est.u.p + est.n.p))
    else log.lik <- sum(log(est.n[cl==0])) + sum(log(est.u.p[cl==1] + est.n.p[cl==1]))
    #cat(log.lik, "\n")
    err <- abs(log.lik.old - log.lik)
    if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
    dec <- ifelse(log.lik.old > log.lik, TRUE, FALSE)
    if(any(pi==0 | pi==1)) err <- 0
  }

  est.u.p <- pi * dunif(x,a,b) 
  est.n.p <- (1-pi) * dnorm(x,mu,sqrt(sigmasq))
  z0 <- est.u.p / (est.n.p + est.u.p)
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
    expr <- (z0 - min(z0)) * sgn.z0

  }
  return(list(expr=expr, a=a, b=b, sigmasq=sigmasq, mu=mu, pi=pi, conv=conv, dec=dec, loc=loc, lik.rec=lik.rec, abnormal=abnormal))
}

