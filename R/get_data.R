### Data Generation ###

# Simulate right-censored processes from the illness-death model with recovery
simulate <- function(n=1000, a12=0.5, a13=1, a21=0.75 ,a23=0.5,
                     p12=1, p13=1, p21=1, p23=1, c.rate=1, b.Z.cens=0.5, 
                     p.Z=0.4, b.Z=0.5, gamma_0=0, gamma_1=1, shape=1,
                     M0=10, M1=50) {
  
  #Simulate the number of patients within each clinic
  m <- rdunif(n=n, min=M0, max=M1)
  
  #Generate cluster id
  cid0 <- unlist(c(sapply(1:n, function(x) rep(x, m[x]))))
  
  #Generate size of the cluster
  M_i0 <- unlist(c(sapply(1:n, function(x) rep(m[x], m[x]))))
  
  #Total number of individuals
  N <- length(cid0)
  
  #Generate individual id
  id0<-1:N
  
  #simulate random effects
  u <- unlist(c(sapply(1:n, function(x) rep(t(rgamma(n = 1, shape=shape)), m[x]))))
  u2 <- unlist(c(sapply(1:n, function(x) rep(t(rgamma(n = 1, shape=shape)), m[x]))))
  u3 <- unlist(c(sapply(1:n, function(x) rep(t(rnorm(n = 1, mean=0, sd=1)), m[x]))))
  
  outdata <- NULL
  for(i in 1:N){
    
    #Simulate covariate that is associated with X(t), right censoring, and missingness 
    Z <- rbinom(n=1, size=1, prob=p.Z)
    
    R = 1
    
    if(R==1){
      
      CC <- rexp(1, rate=c.rate*exp(b.Z.cens*Z)*u2[i])
      
      t12 <- rweibull(1, shape=p12, scale=(1/(a12*exp(b.Z*Z)*u[i])))
      t13 <- rweibull(1, shape=p13, scale=(1/(a13*exp(b.Z*Z)*u[i])))
      stop <- 1*((t13 < t12) | (CC < t12)) #Stop if death or censoring first
      x1 <- min(c(t12, t13)) # Exit time from state 1 (ignoring censoring)
      if (stop==0){
        t1 <- 0
        t2 <- x1
        s1 <- 1
        s2 <- 2
        while (stop==0){
          if(s2[length(s2)]==2){
            t23 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p23, scale=(1/(a23*exp(b.Z*Z)*u[i])))
            t21 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p21, scale=(1/(a21*exp(b.Z*Z)*u[i])))
            
            stop <- 1*((t23 < t21) | (CC < t21)) #Stop if death or censoring first
            x2 <- min(c(t21, t23)) # Exit time from state 2 (ignoring censoring)
            
            t1 <- c(t1, x1)
            t2 <- c(t2, min(x2, CC))
            s1 <- c(s1, 2)
            s2 <- c(s2, ifelse(x2 <= CC, ifelse(t21 <= t23, 1, 3), 2))
          } else {
            t12 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p12, scale=(1/(a12*exp(b.Z*Z)*u[i])))
            t13 <- rtrunc(1, spec="weibull", a=t2[length(t2)], b=Inf,
                          shape=p13, scale=(1/(a13*exp(b.Z*Z)*u[i])))
            
            stop <- 1*((t13 < t12) | (CC < t12)) #Stop if death or censoring first
            x1 <- min(c(t12, t13)) # Exit time from state 1 (ignoring censoring)
            
            t1 <- c(t1, x2)
            t2 <- c(t2, min(x1, CC))
            s1 <- c(s1, 1)
            s2 <- c(s2, ifelse(x1 <= CC, ifelse(t12 <= t13, 2, 3), 1))
          }
        }
      } else {
        t1 <- 0
        t2 <- min(x1, CC)
        s1 <- 1
        s2 <- ifelse(x1 < CC, 3, 1)
      }
    } else {
      t1 <- 0
      t2 <- NA
      s1 <- 0
      s2 <- NA
    }
    
    cid <- rep(cid0[i], length(t1))
    id <- rep(id0[i], length(t1))
    M_i <- rep(M_i0[i], length(t1))
    Z <- rep(Z, length(t1))
    R <- rep(R, length(t1))
    # outdata <- rbind(outdata, data.frame(cid, id, M_i, Tstart = t1, Tstop = t2,
    #                                      from = s1, to = s2, Z, R))
    outdata <- rbind(outdata, data.frame(cid, id, Tstart = t1, Tstop = t2,
                                         from = s1, to = s2, Z, R))
  }
  # t1: start time of the interval
  # t2: end time of the interval
  # s1: state at the start time of the interval
  # s2: state at the end time of the interval
  # if s1==s2 then right censoring occured at t2 at this interval
  outdata <- outdata[order(outdata$cid, outdata$id, outdata$Tstart, outdata$Tstop),]
  return(outdata)
}

reshape_long <- function(dat, tmat) {
  tmat <- matrix(c(NA, 1, 2, 3, NA, 4, NA, NA, NA), nrow = 3, byrow = TRUE)
  new_long <- function(sub_dat, tmat) {
    mat <- lapply(seq_len(ncol(tmat)), function(i) {
      res <- which(!is.na(tmat[i, ]))
      if (length(res)) {
        out <- cbind(from = i, to = res)
        rownames(out) <- NULL
        out
      } else {
        NULL
      }
    })
    mat <- do.call(rbind, mat)
    from <- mat[, 1]
    tmp_dat <- data.frame(cbind(from = mat[, 1], to = mat[, 2],
                                trans = tmat[mat]))
    id <- unique(sub_dat$id)
    is_censor <- sub_dat$from == sub_dat$to
    
    out <- lapply(seq_along(is_censor), function(i) {
      this_censor <- is_censor[i]
      if (this_censor) {
        res <- data.frame(id = id,
                          Tstart = sub_dat$Tstart[i],
                          Tstop = sub_dat$Tstop[i],
                          time = sub_dat$Tstop[i] - sub_dat$Tstart[i],
                          from = sub_dat$from[i])
        res <- merge(res, tmp_dat, by = "from")
        res <- res[, c("id", "from", "to", "Tstart", "Tstop", "time",
                       "trans")]
        res$status <- 0
      } else {
        res <- data.frame(id = id,
                          Tstart = sub_dat$Tstart[i],
                          Tstop = sub_dat$Tstop[i],
                          from = sub_dat$from[i],
                          to = sub_dat$to[i],
                          time = sub_dat$Tstop[i] - sub_dat$Tstart[i],
                          status = 1)
        res <- merge(res, tmp_dat, by = c("from", "to"))
        res <- res[, c("id", "from", "to", "trans", "Tstart",
                       "Tstop", "time", "status")]
      }
      res
    })
    do.call(rbind, out)
  }
  res <- by(data = dat, dat$id, new_long, tmat = tmat)
  res <- do.call(rbind, res)
  class(res) <- c( "data.frame")
  merge_cid <- data.frame( id = dat$id, cid = dat$cid)
  res <- merge(merge_cid, res, by = "id")
  return(res)
}


#######################################
#### Add group information#############
#######################################

add_group <- function(data){
  group <- by(data, data$cid, function(sub_dat){
    l <- length(unique(sub_dat$id))
    rbinom(l, 1, 0.5) 
  })
  g <- do.call(c, group)
  ## Add group information
  res <- do.call(rbind, by(data, data$id, function(sub_dat){
    sub_dat$group <- rep(g[sub_dat$id[1]], nrow(sub_dat))
    return(sub_dat)
  }))
  return(res)
}

trans <- function(nstate, state_names, from, to) {
  if (missing(nstate) && missing(state_names))
    stop("One of 'nstate' and 'state_names' has to be specified.")
  if (missing(state_names)) {
    state_names <- as.character(seq_len(nstate))
  } else {
    state_names <- unique(state_names)
    nstate <- length(state_names)
  }
  if (length(from) != length(to))
    stop("The length of 'from' and 'to' must be the same.")
  if (is.character(from)) {
    from <- match(from, state_names)
  } else {
    from <- as.integer(from)
  }
  if (is.character(to)) {
    to <- match(to, state_names)
  } else {
    to <- as.integer(to)
  }
  mat <- matrix(FALSE, ncol = nstate, nrow = nstate)
  dimnames(mat) <- list(state_names, state_names)
  mat[cbind(from, to)] <- TRUE
  mat
}


get_data <- function(n, M0 = 10, M1 = 20, tmat){
  dat <- simulate(n = n, M0 = M0, M1 = M1)
  dat <- reshape_long(dat, tmat = tmat)
  dat <- add_group(dat)
  return(dat)
}
