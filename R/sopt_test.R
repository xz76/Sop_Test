

sopt_test <- function(data, tmat, cid = "cid", id = "id", group = "group", 
                      j = 2, B = 1000, ipw = 0, tau = NULL, method = "linear"){
  
  Pt <- list()
  tms <- list()
  S <- as.numeric(sort(unique(c(names(table(data$from)), names(table(data$to))))))
  T_c <- as.numeric(sort(unique(names(table(data$from)))))
  
  if(is.null(data[[group]])){
    stop("group information is needed")
  }

  if(nrow(tmat)!=ncol(tmat)){
    stop("tmat is not a square matrix")
  }
  if(min(S %in% 1:nrow(tmat))==0){
    stop("states should be numbered as 1 to total number of states")
  }
  if(min(S[rowSums(tmat)==0] %in% T_c)){
    stop("Absorbing states according to tmat appeared in the from variable")
  }
  if(!(method %in% c("linear", "KS"))){
    stop("Possible choices for method 'linear' and 'KS'")
  }
  if(B<=0){
    stop("Tests cannot be performed based on <=0 bootstrap samples")
  } else if (method=="KS" & B<1000){
    warning("It is recommended to use at least 1000 bootstrap samples for the Kolmogorov-Smirnov-type test")
  } else if (method=="linear" & B<100){
    warning("It is recommended to use at least 100 bootstrap samples for the linear test")
  }
  n <- length(unique(data[,cid]))
  groups <- unique(data[,group])
  groups <- sort(groups[!is.na(groups)])
  if(length(groups)!=2){
    stop("Number of groups != 2")
  }
  
  for(g in groups){
      group_data <- data[data[,group]==g,]
      P0 <- sop_t(data = group_data, tau=NULL,
                  ipw=0, tmat = tmat, times = NULL)
    rep1 <- function(x){
      x[c(1, seq_along(x))]
    }
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[1]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
      
    } else {
      Pt[[2]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[2]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
    }
  }
  
  tms <- sort(unique(c(tms[[1]], tms[[2]])))
  
  if(nrow(data[data$from==j,])>0){
    tS <- sort(unique(c(data[data$to==j,"from"],j)))
  } else {
    tS <- sort(unique(data[data$to==j,"from"]))
  }
  
  EY <- NULL
  EY.t <- function(t,dt){
    sum(dt$Tstart<t & dt$Tstop>=t)/n
  }
  for(i in tS){
    for(g in groups){
      dat <- unique(data[data$from==i & data[,group]==g,
                         c(id, "from", "Tstart", "Tstop")])
      if(is.null(EY)){
        EY <- sapply(tms, EY.t, dt=dat)
      } else {
        EY <- cbind(EY, sapply(tms, EY.t, dt=dat))
      }
    }
  }
  
  Wt <- rowProds(EY)/rowSums(EY)
  
  tms <- tms[!is.na(Wt)]
  if(length(tms)==0 | max(Wt, na.rm=TRUE)==0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]

  
  estimator <- function(data,cov,tau,ipw, tmat, Wt, method){
    if (cov != 0 && cov != 1) {
      stop("The 'cov' has to be 0 or 1.")
    }
    
    P1 <- sop_t(data[cov==1,], tau=tau, ipw=ipw, tmat = tmat)
    P0 <- sop_t(data[cov==0,], tau=tau,  ipw=ipw, tmat = tmat)
    
    p1 <- P1[,c("t",paste("p",j,sep=""))]
    p0 <- P0[,c("t",paste("p",j,sep=""))]
    
    elemenTstart<-function(t){
      max(1:length(p1$t)*(p1$t<=t))
    }
    elementfrom<-sapply(tms,elemenTstart)
    elementfrom<-(elementfrom==0)+(elementfrom>0)*elementfrom
    
    element0<-function(t){
      max(1:length(p0$t)*(p0$t<=t))
    }
    elements0<-sapply(tms,element0)
    elements0<-(elements0==0)+(elements0>0)*elements0
    
    D_t <- p1[elementfrom,paste("p",j,sep="")]-p0[elements0,paste("p",j,sep="")]
    tau0 <- max(data$Tstop)
    tmax <- min(max(p1[p1[,2]>0, "t"]),
                max(p0[p0[,2] >0, "t"]))
    tau <- min(tau0, tmax)
    
    if (method == "linear"){
      dm<-diff(c(tms,tau),lag=1)
      res <- sum(Wt * D_t * dm)
    } else if(method=="KS"){
      D_hat <- (Pt[[1]](tms) - Pt[[2]](tms))
      res <- max(abs(Wt * (D_t - D_hat)))
    }
    return(res)
  }
  dauc_boot <- function (data, cov, tau =NULL, ipw, tmat,
                         Wt, B,  verbose = 0, id,  method)
  {
    ids <- unique(data[[id]])
    n <- length(ids)
    res <- matrix(NA, length(Z0), B)
    
    for (b in 1:B) {
      if (verbose > 0) {
        cat("\nBootstrap replication", b, "\n")
        flush.console()
      }
      bootdata <- NULL
      bids <- sample(ids, replace = TRUE)
      bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
      bootdata <- data[bidxs, ]
      if (verbose > 0) {
        print(date())
        print(events(bootdata))
        cat("applying theta ...")
      }
      if (all(bootdata$from == bootdata$to)){
        next
      }
      thstar <-  try(estimator( data = bootdata, cov = bootdata$group, tau = tau,
                                ipw = ipw, tmat = tmat, Wt = Wt, method = method))
      if (class(thstar)!="try-error"){
        res[, b] <- thstar
      }else{
        res[, b] <- NA
      }
      
    }
    if (verbose)
      cat("\n")
    return(res)
  }
  Z0 <-  estimator( data = data, cov = data[, group], tau = NULL, 
                    ipw = 0, tmat = tmat, Wt = Wt,  method = method)
  if (method == "linear") {

    dauc_res <- dauc_boot( data = data, cov = data$group,
                           ipw = ipw, tmat = tmat, Wt = Wt, B = B, id = cid, method = method)
    dauc_res <- as.numeric(na.omit(dauc_res))
    dauc_sd <- sd(dauc_res)
    
    T_linear <- Z0/dauc_sd
    pval <- 2*pnorm(abs(T_linear), lower.tail = FALSE)
  } else {
    D_hat <- (Pt[[1]](tms) - Pt[[2]](tms))
    KS0 <- max(abs(Wt*D_hat)* sqrt(n))
    browser()
    KS.b <- dauc_boot(data = data, cov = data$group,
                      ipw = ipw, tmat = tmat, Wt = Wt, B = B, id = cid, method = method)
    KS.boot <- apply(abs(KS.b),2,max)
    pval <- mean(KS.boot >= KS0, na.rm = TRUE)
  }
  
  names(pval) <- paste0("p-value at State", j)
  return(pval)
}




