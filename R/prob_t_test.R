
tpr.test <- function(data, tmat, cid = "cid", id = "id", group = "group", 
                      j = 2, B = 1000, s = NULL, h = NULL, tau = NULL,
                     method = "linear"){
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
  
  n <- length(unique(data[[cid]]))
  groups <- unique(data[[group]])
  if(length(groups)!=2){
    stop("Number of groups != 2")
  }
  groups <- c(0, 1)
  data[[group]] <- factor(data[[group]], labels = c(0, 1))
  data[[group]] <- as.integer(data[[group]]) - 1
  
  lab <- ifelse(is.null(h), paste0("p", j), paste0("p", h, j))
  for(g in groups){
    group_data <- data[data[[group]] == g, ]
    P0 <- tpr.t(data = group_data, s = s, h = h, j = j, tau = tau,
                tmat = tmat, times = NULL, CI = FALSE)
  
    rep1 <- function(x){
      x[c(1, seq_along(x))]
    }
    
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$t, rep1(P0[, lab]))
      tms[[1]] <- P0$t[diff(rep1(P0[, lab]), lag = 1) != 0]
      
    } else {
      Pt[[2]] <- stepfun(P0$t, rep1(P0[, lab]))
      tms[[2]] <- P0$t[diff(rep1(P0[, lab]), lag = 1) != 0]
    }
  }
  
  tms <- sort(unique(c(tms[[1]], tms[[2]])))
  
  if(nrow(data[data$from==j,])>0){
    tS <- sort(unique(c(data[data$to == j, "from"], j)))
  } else {
    tS <- sort(unique(data[data$to == j, "from"]))
  }
  
  EY <- NULL
  EY.t <- function(t,dt){
    sum(dt$Tstart < t & dt$Tstop >= t)/n
  }
  for(i in tS){
    for(g in groups){
      dat <- unique(data[data$from == i & data[[group]] == g,
                         c(id, "from", "Tstart", "Tstop")])
      if(is.null(EY)){
        EY <- sapply(tms, EY.t, dt = dat)
      } else {
        EY <- cbind(EY, sapply(tms, EY.t, dt = dat))
      }
    }
  }
  
  Wt <- rowProds(EY)/rowSums(EY)
  
  tms <- tms[!is.na(Wt)]
  if(length(tms) == 0 | max(Wt, na.rm = TRUE) == 0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]

  estimator <- function(data, cov, tau, s, h, j, tmat, Wt, D_hat = NULL){
    if (cov != 0 && cov != 1) {
      stop("The 'cov' has to be 0 or 1.")
    }
    P1 <- tpr.t(data[cov == 1, ], tau=tau, s=s, h=h, j=j, 
                tmat = tmat, times=NULL, CI = FALSE)
    P0 <- tpr.t(data[cov == 0, ], tau=tau, s=s, h=h, j=j, 
                tmat = tmat, times=NULL, CI = FALSE)
    p1 <- P1[, c("t", lab)]
    p0 <- P0[, c("t", lab)]
    
    elemenTstart < -function(t){
      max(1:length(p1$t)*(p1$t <= t))
    }
    elementfrom <- sapply(tms, elemenTstart)
    elementfrom <- (elementfrom == 0) + (elementfrom > 0) * elementfrom
    
    element0<-function(t){
      max(1 : length(p0$t) * (p0$t <= t))
    }
    elements0 <- sapply(tms, element0)
    elements0 <- (elements0 == 0) + (elements0 > 0) * elements0
    
    D_t <- p1[elementfrom, lab] - p0[elements0, lab]
    if (method == "linear"){
      tau0 <- max(data$Tstop)
      tmax <- min(max(p1[p1[ , 2] > 0, "t"]),
                  max(p0[p0[ , 2] > 0, "t"]))
      tau <- min(tau0, tmax)
      dm<-diff(c(tms, tau), lag = 1)
      res <- sum(Wt * D_t * dm)
    } else if(method == "KS"){
      res <- max(abs(Wt * (D_t - D_hat)))
    }
    return(res)
  }
  dauc_boot <- function (data, cov, tau =NULL,  tmat,
                         Wt, B, id, D_hat=NULL, h = h, j = j, s = s)
  {
    ids <- unique(data[[id]])
    n <- length(ids)
    res <- matrix(NA, 1, B)
    for (b in 1:B) {
      bootdata <- NULL
      bids <- sample(ids, replace = TRUE)
      bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
      bootdata <- data[bidxs, ]
      if (all(bootdata$from == bootdata$to)){
        next
      }
      thstar <-  try(estimator( data = bootdata, cov = cov, tau = tau,
                               tmat = tmat, Wt = Wt, D_hat = D_hat, h = h, 
                               j = j, s = s))
      if (class(thstar)!="try-error"){
        res[ , b] <- thstar
      }else{
        res[ , b] <- NA
      }
    }
    return(res)
  }

  if (method == "linear") {
    Z0 <-  estimator( data = data, cov = data[[group]], tau = NULL, 
                      tmat = tmat, Wt = Wt, h = h, j = j, s = s)
    dauc_res <- dauc_boot( data = data, cov =  data[[group]],
                          tmat = tmat, Wt = Wt, B = B, id = cid, 
                          h = h, j = j, s = s)
    dauc_res <- as.numeric(dauc_res)
    dauc_sd <- sd(dauc_res, na.rm = TRUE)
    T_linear <- Z0 / dauc_sd
    pval <- 2 * pnorm(abs(T_linear), lower.tail = FALSE)
  } else {
    D_hat <- (Pt[[1]](tms) - Pt[[2]](tms))
    KS0 <- max(abs(Wt * D_hat))
    KS.b <- dauc_boot(data = data, cov =  data[[group]],
                     tmat = tmat, Wt = Wt, B = B, id = cid,
                     D_hat = D_hat, h = h, j = j, s = s)
    pval <- mean(KS.b >= KS0, na.rm = TRUE)
  }
  ifelse(is.null(s), names(pval) <- paste0("p-value at State", j) ,
         names(pval) <- paste0("p-value at P",h, j))
  return(pval)
}




