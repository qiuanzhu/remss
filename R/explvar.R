explvar <- function(par, datalines, ioparms) {
    
    n <- ioparms[1]
    nk <- ioparms[2]
    np <- ioparms[3]
    
    zeit <- datalines[, 1]
    exbeta <- rep(0, ioparms[1])
    istat <- datalines[, 2]
    inde <- rep(0, ioparms[1])
    zeitp <- datalines[, nk + 3]
    s0 <- datalines[, nk + 4]
    skm <- datalines[, nk + 5]
    bkm <- datalines[, nk + 6]
    dup <- rep(0, ioparms[3])
    cov <- datalines[, 3:(nk + 2)]
    exbeta <- rep(0, n)
    iz <- 0
    
    for (j in 1:np) {
        dup[j] <- 0
        
        if (iz == n) {
            break
        }
        iz <- iz + 1
        
        while (istat[iz] == 0 & iz != n) {
            inde[iz] <- j
            iz <- iz + 1
        }
        if (iz == n) {
            break
        }
        
        while (zeit[iz] == zeitp[j] & iz != n) {
            dup[j] <- dup[j] + 1
            inde[iz] <- j
            iz <- iz + 1
            while (istat[iz] == 0 & iz != n) {
                inde[iz] <- j
                iz <- iz + 1
            }
        }
        if (zeit[iz] > zeitp[j]) {
            iz <- iz - 1
        }
    }
    
    d1u_s <- 0
    d1c_s <- 0
    wei_s <- 0
    
    for (i in 1:n) {
        xbeta <- 0
        for (j in 1:nk) {
            xbeta <- xbeta + cov[i, j] * par[j]
        }
        exbeta[i] <- exp(xbeta)
    }
    
    
    for (j in 1:np) {
        resc <- 0
        resk <- 0
        for (i in 1:n) {
            cox <- s0[j]^exbeta[i]
            
            if (zeit[i] > zeitp[j]) {
                resc <- resc + (1 - cox)
                resk <- resk + (1 - skm[j])
            } else if (zeit[i] <= zeitp[j] & istat[i] != 0) {
                resc <- resc + cox
                resk <- resk + skm[j]
            } else {
                ife <- inde[i]
                coxc <- s0[ife]^exbeta[i]
                if (coxc > 1e-06) {
                  resc <- resc + (1 - cox) * cox/coxc + cox * (1 - cox/coxc)
                }
                if (skm[ife] > 1e-06) {
                  resk <- resk + (1 - skm[j]) * skm[j]/skm[ife] + skm[j] * (1 - skm[j]/skm[ife])
                }
            }
            
        }
        d1u_s <- d1u_s + (resk/n) * dup[j]/bkm[j]
        d1c_s <- d1c_s + (resc/n) * dup[j]/bkm[j]
        wei_s <- wei_s + dup[j]/bkm[j]
        
    }
    
    
    
    ioparms[4] <- d1u_s/wei_s
    ioparms[5] <- d1c_s/wei_s
    ioparms[6] <- (ioparms[4] - ioparms[5])/ioparms[4]
    
    return(ioparms[6])
}
