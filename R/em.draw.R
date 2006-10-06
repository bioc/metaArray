"em.draw" <-
function(vec, cl=NULL, threshold=0.0001) {
  par(mfrow=c(2,2))
  raw <- as.numeric(vec) 
  if(is.null(cl)) cl <- rep(0,length(vec))
  fit.tmp <- fit.em(raw, cl, threshold)
  trans <- fit.tmp$expr
  lik <- fit.tmp$lik.rec
  ord <- order(raw)
  raw <- raw[ord]; trans <- trans[ord]
  plot(raw, trans, main="estimated poe", xlab="raw expression", ylab="probability", ylim=c(-1.1,1.1), pch=3,col=2)
  lines(raw, trans, lty=3, col=2)
  abline(0,0, lty=2, col=4)
  mu <-fit.tmp$mu;sigmasq<-fit.tmp$sigmasq
  pi<-fit.tmp$pi
  a <-fit.tmp$a; b<-fit.tmp$b
  xvals <- seq(a,b,by=0.01)
  deu <- dunif(xvals,a,b)
  den <- dnorm(xvals,mu,sigmasq)
  dem <- pi * deu + (1-pi) * den
  plot(xvals, dem, type="n", main="mixture density", xlab="expression", ylab="density", ylim=c(0, max(c(den,dem,deu))*1.1))
  lines(xvals, dem, lty=2, col=2)
  lines(xvals, deu, lty=3, col=4)
  lines(xvals, den, lty=3, col=4)
  hist(raw, main="raw expression")
  lik.max <- max(lik) 
  lik.min <- min(lik)
  lik.len <- length(lik)
  plot(1:lik.len, lik, xlim=c(0,lik.len+1), ylim=c(lik.min,lik.max), main="log likelihood", xlab="iteration", ylab="likelihood", type="n")
  lines(1:lik.len, lik, lty=1, col=2)
}

