
"poe.mcmc" <-
function (AA, NN = NULL, id = NULL, M = 2000, kap.min=3.0,
    logdata=FALSE,     stepsize = 0.5, 
    centersample = FALSE, centergene = FALSE, generatestarts = TRUE, start.method = 1, 
    startobject = R0, collapse.to.two = FALSE, burnin = 200,
    collapse.window=50,converge.threshold=0.01,
    PR = list(alpha.mm = 0, alpha.sd = 100, mu.mm = 0, mu.sd = 100, 
        pipos.mm = 0, pipos.sd = 100, pineg.mm = 0, pineg.sd = 100, 
        kap.pri.rate = 1, tausqinv.aa = 1, tausqinv.bb = 0.1)) 
{

    TT <- ncol(AA)
    GG <- nrow(AA)
    tmp.poemat <- matrix(0,GG,TT)
    if (!is.null(NN) & length(NN) != TT) {
        stop("Length of NN should be the same as ncol(AA)")
    }
    if(is.na(sum(as.vector(AA))))
      stop("NA values not allowed in AA matrix. Remove rows with NA values and retry.")
    if (logdata == TRUE) 
        AA <- log(AA)
    if (is.null(id) == TRUE) 
        id <- as.character(1:GG)
    if (is.null(NN)) 
        NN <- rep(0, TT)
    if (centersample) {
        sammed <- apply(AA, 2, median)
        AA <- sweep(AA, 2, sammed)
    }
    if (centergene) {
        genmed <- apply(AA, 1, median)
        AA <- sweep(AA, 1, genmed)
    }
    if (generatestarts == TRUE) {
        alpha0 <- apply(AA, 2, mean)
        mats <- t(apply(AA, 1, sort))
        mug0 <- sigmag0 <- kappaposg0 <- kappanegg0 <- rep(NA, 
            GG)
        piposg0 <- pinegg0 <- rep(NA, GG)
        if (start.method == 1) {
            for (gg in 1:GG) {
                mug0[gg] <- mean(mats[gg, -c(1, 2, TT - 1, TT)])
                pinegg0[gg] <- piposg0[gg] <- 2/TT
                sigmag0[gg] <- sqrt(var(mats[gg, -c(1, 2, TT - 
                  1, TT)]))
                kappaposg0[gg] <- max(kap.min * sigmag0[gg], 
                  mug0[gg] - mats[gg, 1])
                kappanegg0[gg] <- max(kap.min * sigmag0[gg], 
                  mats[gg, TT] - mug0[gg])
            }
        }
        if (start.method == 2) {
            for (gg in 1:GG) {
                clu <- kmeans(mats[gg, ], 2)
                if (clu$size[1] == clu$size[2]) {
                  big <- 1
                  sma <- 2
                }
                else {
                  big <- (1:2)[clu$size == max(clu$size)]
                  sma <- (1:2)[clu$size == min(clu$size)]
                }
                mug0[gg] <- clu$centers[big]
                sigmag0[gg] <- sqrt(var(mats[gg, clu$cluster == 
                  big]))
                if (mug0[gg] > 0) {
                  pinegg0[gg] <- clu$size[sma]/TT
                  piposg0[gg] <- 1/TT
                }
                else {
                  piposg0[gg] <- clu$size[sma]/TT
                  pinegg0[gg] <- 1/TT
                }
                kappaposg0[gg] <- max(max(mats[gg, ]) - mug0[gg], 
                  kap.min * sigmag0[gg])
                kappanegg0[gg] <- max(mug0[gg] - min(mats[gg, 
                  ]), kap.min * sigmag0[gg])
            }
        }
        if (start.method == 3) {
            for (gg in 1:GG) {
                clu <- kmeans(mats[gg, ], 3)
                posc <- (1:3)[clu$centers == max(clu$centers)]
                negc <- (1:3)[clu$centers == min(clu$centers)]
                midc <- (1:3)[clu$centers == median(clu$centers)]
                mug0[gg] <- clu$centers[midc]
                sigmag0[gg] <- sqrt(var(mats[gg, clu$cluster == 
                  midc]))
                pinegg0[gg] <- clu$size[negc]/TT
                piposg0[gg] <- clu$size[posc]/TT
                kappaposg0[gg] <- max(max(mats[gg, ]) - mug0[gg], 
                  kap.min * sigmag0[gg])
                kappanegg0[gg] <- max(mug0[gg] - min(mats[gg, 
                  ]), kap.min * sigmag0[gg])
            }
        }
        if (start.method == 4) {
            if (sum(NN) < 3) 
                stop("start.method==4 requires at least 3 normal samples")
            for (gg in 1:GG) {
                mug0[gg] <- mean(AA[gg, NN == 1])
                sigmag0[gg] <- sqrt(var(AA[gg, NN == 1]))
                clu <- kmeans(AA[gg, ], 2)
                norm <- (1:2)[abs(clu$centers - mug0[gg]) == 
                  min(abs(clu$centers - mug0[gg]))]
                expr <- 1
                if (norm == 1) 
                  expr <- 2
                if (clu$centers[norm] > clu$centers[expr]) {
                  pinegg0[gg] <- min(clu$size[expr]/TT, 1/2)
                  piposg0[gg] <- 1/TT
                }
                else {
                  piposg0[gg] <- min(clu$size[expr]/TT, 1/2)
                  pinegg0[gg] <- 1/TT
                }
                kappaposg0[gg] <- max(max(AA[gg, ]) - mug0[gg], 
                  kap.min * sigmag0[gg])
                kappanegg0[gg] <- max(mug0[gg] - min(AA[gg, ]), 
                  kap.min * sigmag0[gg])
            }
        }
        sigmagsqinv0 <- 1/sigmag0^2
        gamma0 <- (mean(sigmagsqinv0))^2/var(sigmagsqinv0)
        lambda0 <- mean(sigmagsqinv0)/var(sigmagsqinv0)
        tmpmat <- matrix(0,GG,TT)

        PP <- list(alpha = alpha0, mug = mug0, kappaposg = kappaposg0, 
            kappanegg = kappanegg0, sigmag = sigmag0, piposg = piposg0, 
            pinegg = pinegg0, mu = mean(mug0), tausqinv = 1/var(mug0), 
            gamma = gamma0, lambda = lambda0, pil.pos.mean = mean(logit(piposg0)), 
            pil.neg.mean = mean(logit(pinegg0)), pil.pos.prec = mean(piposg0) * 
                (1 - mean(piposg0)), pil.neg.prec = mean(pinegg0) * 
                (1 - mean(pinegg0)), kap.pos.rate = 1/mean(kappaposg0), 
            kap.neg.rate = 1/mean(kappanegg0), poe=tmp.poemat, phat.pos=tmpmat, 
            phat.neg=tmpmat, accept = 0)
    }
    if (generatestarts == FALSE) {
        mmm <- length(startobject$gamma)
        PP <- list(alpha = startobject$alpha[mmm, ], mug = startobject$mug[mmm, 
            ], kappaposg = startobject$kappaposg[mmm, ], kappanegg = startobject$kappanegg[mmm, 
            ], sigmag = startobject$sigmag[mmm, ], piposg = startobject$piposg[mmm, 
            ], pinegg = startobject$pinegg[mmm, ], mu = startobject$mu[mmm], 
            tausqinv = startobject$tausqinv[mmm], gamma = startobject$gamma[mmm], 
            lambda = startobject$lambda[mmm], pil.pos.mean = startobject$pil.pos.mean[mmm], 
            pil.pos.prec = startobject$pil.pos.prec[mmm], pil.neg.mean = startobject$pil.neg.mean[mmm], 
            pil.neg.prec = startobject$pil.neg.prec[mmm], kap.pos.rate = startobject$kap.pos.rate[mmm], 
            kap.neg.rate = startobject$kap.neg.rate[mmm], poe=tmp.poemat, 
            phat.pos=tmpmat, phat.neg=tmpmat, accept = 0)
    }
  if (collapse.to.two) {
    cat("Collapsing to Two: ")
    cw <- collapse.window
    mm <- 0
      pidif <- 1
      pim.pos <- pim.neg <- 100
      while ( pidif>converge.threshold) {
        mm <- mm+1
         new.PR <- unlist(PR)
         new.PP <- unlist(PP)
         avgpos <- rep(0,GG*(3*TT+6)+TT+11)
         res <- .C("poe_fit_2",as.matrix(AA),as.integer(NN),
            as.double(new.PR), as.double(new.PP),
            as.integer(GG), as.integer(TT),as.integer(1),
            avgpos=as.double(avgpos), PACKAGE="metaArray")$avgpos
         PP$alpha <- res[1:TT]
         PP$mug <- res[(TT+1):(TT+GG)]
         PP$kappaposg <- res[(TT+GG+1):(TT+2*GG)]
         PP$kappanegg <- res[(TT+2*GG+1):(TT+3*GG)]
         PP$sigmag <- res[(TT+3*GG+1):(TT+4*GG)]
         PP$piposg <- res[(TT+4*GG+1):(TT+5*GG)]
         PP$pinegg <- res[(TT+5*GG+1):(TT+6*GG)]
         PP$mu <- res[TT+6*GG+1]
         PP$tausqinv <- res[TT+6*GG+2]
         PP$gamma <- res[TT+6*GG+3]
         PP$lambda <- res[TT+6*GG+4]
         PP$pil.pos.mean <- res[TT+6*GG+5]
         PP$pil.neg.mean <- res[TT+6*GG+6]
         PP$pil.pos.prec <- res[TT+6*GG+7]
         PP$pil.neg.prec <- res[TT+6*GG+8]
         PP$kap.pos.rate <- res[TT+6*GG+9]
         PP$kap.neg.rate <- res[TT+6*GG+10]
         PP$poe <- matrix(res[(TT+6*GG+11):(GG*(TT+6)+TT+10)], nrow=GG)
         PP$phat.pos <- matrix(res[(GG*(TT+6)+TT+11):(GG*(2*TT+6)+TT+10)], nrow=GG)
         PP$phat.neg <- matrix(res[(GG*(2*TT+6)+TT+11):(GG*(3*TT+6)+TT+10)], nrow=GG)
         PP$accept <- res[GG*(3*TT+6)+TT+11]
         gc();
         pim.pos <- c(pim.pos, mean(PP$piposg))
         pim.neg <- c(pim.neg, mean(PP$pinegg))
        
        if(mm>(2*cw)) {
          pidif <- 0.5*max(
                  abs(mean(pim.neg[(mm-cw-1  ):(mm   )])-
                      mean(pim.neg[(mm-2*cw-1):(mm-cw)]))/
                     (mean(pim.neg[(mm-cw-1  ):(mm   )])+
                      mean(pim.neg[(mm-2*cw-1):(mm-cw)])),
                  abs(mean(pim.pos[(mm-cw-1  ):(mm   )])-
                      mean(pim.pos[(mm-2*cw-1):(mm-cw)]))/
                     (mean(pim.pos[(mm-cw-1  ):(mm   )])+
                      mean(pim.pos[(mm-2*cw-1):(mm-cw)]))) 
         }
     }
     for (gg in 1:GG) {
           ppplus <- mean(PP$phat.pos[gg, ])
           ppminu <- mean(PP$phat.pos[gg, ]-PP$poe[gg,])
           pp0 <- 1 - ppminu - ppplus
           if (ppplus < 2/TT & abs(ppminu - pp0) < 6/TT) {
               PP$phat.pos[gg, ] <- abs(PP$poe[gg, ])
               PP$poe[gg,] <-abs(PP$poe[gg,])
               ee <- round(PP$poe[gg,])
               PP$mug[gg] <- mean(AA[gg,ee==0],na.rm=TRUE)
               PP$piposg[gg] <- mean(ifelse(ee==1,1,0))
               PP$pinegg[gg] <- mean(ifelse(ee==-1,1,0))
               PP$sigmag[gg] <- sqrt(var(AA[gg,ee==0],na.rm=TRUE))
               PP$kappanegg[gg] <- kap.min*PP$sigmag[gg]
               PP$kappaposg[gg] <- kap.min*PP$sigmag[gg]
          }
     }  
     cat("Finished\n")
  }  
  cat("Gibbs Sampler\n")
  new.PR <- unlist(PR)
  new.PP <- unlist(PP)
  avgpos <- rep(0,GG*(3*TT+6)+TT+11)
  res <- .C("poe_fit",as.matrix(AA),as.integer(NN),
       as.double(new.PR),as.double(new.PP),as.integer(GG),
       as.integer(TT),as.integer(M),as.integer(burnin),
       avgpos=as.double(avgpos), PACKAGE="metaArray")$avgpos
  alpha <- res[1:TT]
  mug <- res[(TT+1):(TT+GG)]
  kappaposg <- res[(TT+GG+1):(TT+2*GG)]
  kappanegg <- res[(TT+2*GG+1):(TT+3*GG)]
  sigmag <- res[(TT+3*GG+1):(TT+4*GG)]
  piposg <- res[(TT+4*GG+1):(TT+5*GG)]
  pinegg <- res[(TT+5*GG+1):(TT+6*GG)]
  mu <- res[TT+6*GG+1]
  tausqinv <- res[TT+6*GG+2]
  gamma <- res[TT+6*GG+3]
  lambda <- res[TT+6*GG+4]
  pil.pos.mean <- res[TT+6*GG+5]
  pil.neg.mean <- res[TT+6*GG+6]
  pil.pos.prec <- res[TT+6*GG+7]
  pil.neg.prec <- res[TT+6*GG+8]
  kap.pos.rate <- res[TT+6*GG+9]
  kap.neg.rate <- res[TT+6*GG+10]
  poe <- matrix(res[(TT+6*GG+11):(GG*(TT+6)+TT+10)], nrow=GG)
  poe <- as.data.frame(poe)
  rownames(poe) <- rownames(AA)
  colnames(poe) <- colnames(AA)
  accept <- res[GG*(3*TT+6)+TT+11]
  gc();
  foo <- function(org, fit, inter = 0.01) {
    ord <- order(org)
    org.t <- org[ord]
    fit.t <- fit[ord]
    l <- length(org)
    min.fit <- min(fit.t)
    i <- 1
    while(fit.t[i]-min.fit > inter) {
      fit.t[i] <- min.fit
      i <- i+1
    }
    i <- length(org)
    max.fit <- max(fit.t)
    while(max.fit - fit.t[i] > inter) { 
      fit.t[i] <- max.fit
      i <- i-1
    }
    fit2 <- fit
    fit2[ord] <- fit.t
    fit2
  }
  for(i in 1:GG) {
    poe[i,] <- foo(as.numeric(AA[i,]), as.numeric(poe[i,]))
  }
  return(list(alpha = alpha, mug = mug, kappaposg = kappaposg, 
        kappanegg = kappanegg, sigmag = sigmag, piposg = piposg, 
        pinegg = pinegg, mu = mu, tausqinv = tausqinv, gamma = gamma, 
        lambda = lambda, pil.pos.mean = pil.pos.mean, pil.pos.prec = pil.pos.prec, 
        pil.neg.mean = pil.neg.mean, pil.neg.prec = pil.neg.prec, 
        kap.pos.rate = kap.pos.rate, kap.neg.rate = kap.neg.rate,
        poe=poe, accept = accept))
}

