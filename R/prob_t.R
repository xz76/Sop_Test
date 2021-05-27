prob.est <- function(data, s=NULL, h=NULL, j, tau = NULL, tmat,
                  times = NULL, weights = NULL, CI = TRUE){
  S <- as.numeric(sort(unique(c(data$from, data$to))))
  T_c <- as.numeric(sort(unique(data$from)))
  if(nrow(tmat) != ncol(tmat)){
    stop("tmat is not a square matrix")
  }
  if(min(S %in% 1 : nrow(tmat)) == 0){
    stop("states should be numbered as 1 to total number of states")
  }
  if(min(S[rowSums(tmat) == 0] %in% T_c)){
    stop("Absorbing states according to tmat appeared in the from variable")
  }
  if(!(j %in% S)){
    stop("State j not in the state space of the process")
  }
  if(is.null(h) & !is.null(s)){
    warning("The starting state h is null. State occupation probabilities are
            estimated which ignore the argument s")
  }
  if(is.null(tau)){
    tau <- max(data$Tstop)
  }
  
  ## unique time points
  tt <-sort(unique(data[data$to != data$from,"Tstop"]))
  tt <- tt[tt<=tau]
  if(!is.null(s)){
    tt <- tt[tt>s]
  }
  ## prob_n function to estimate probability at state j   
  prob_n <- function(data){
  CTI <- list()
  counter <- 1
  for(h0 in T_c){
    for(j0 in S[tmat[h0, ]]){
      data_h <- data[data$from == h0, ]
      data_h$delta <- 1*(data_h$to == j0)
      if(is.null(weights)){
        fit <- coxph(Surv(Tstart, Tstop, delta, type = "counting") ~ 1, 
                     data = data_h, control = coxph.control(timefix = FALSE))
      } else {
        fit <- coxph(Surv(Tstart, Tstop, delta, type = "counting") ~ 1,
                     weight = weights, data = data_h,
                     control = coxph.control(timefix = FALSE))
      }
      A <- basehaz(fit, centered=FALSE)
      A_t<-stepfun(A$time, c(0, A$hazard))
      CTI[[counter]] <- A_t
      if(counter==1){
        pointer <- c(h0, j0, counter)
      } else {
        pointer <- rbind(pointer, c(h0, j0, counter))
      }
      counter <- counter + 1
    }
  }
  dA <- sapply(seq_along(CTI), function(i) {
    diff(c(0, CTI[[i]](tt)), lag = 1)
  })
  ## Corresponding to pointer.  
  ttrans <- t(tmat)
  mat0 <- matrix(0, nrow = nrow(tmat), ncol = ncol(tmat))
  mat_list <- lapply(seq_len(nrow(dA)), function(i) {
    out <- mat0
    out[which(ttrans)] <- dA[i, ]
    out <- t(out)
    ## compute diagonal elements
    diag(out) <- - rowSums(out)
    out + diag(nrow(tmat))
  })
  
  ## Calculate the product integral estimator
  P_n <-Reduce("%*%",  mat_list, accumulate = TRUE)
  
    if(is.null(h)){
      p_0 <- sapply(S, function(i){
        sum(data[data$Tstart==0, "from"] == i) / nrow(data[data$Tstart == 0, ])
      })
      p_n <- matrix(NA, nrow = length(tt), ncol = nrow(tmat))
      for(i in 1 : length(tt)){
        p_n[i,] <- p_0 %*% P_n[[i]]
      }
      p_n <- rbind(c(0, p_0), cbind(tt, p_n))
      p_n <- as.data.frame(p_n)
      colnames(p_n) <- c("t", paste0("p", S))
      rownames(p_n) <- 1:(length(tt)+1)
      p_n <- p_n[,c("t", paste0("p", j))]
    } else {
      p_n <- unname(sapply(P_n, function(x) x[h, j]))
      p_n <- as.data.frame(cbind(tt, p_n))
      colnames(p_n) <- c("t",  paste0("p", h, j))
    } 
    if (!is.null(times)){
      p_nt <- sapply(times, function(t){tail(p_n[which(p_n$t < t), ], 1)})
      colnames(p_nt) <- times
      res <- t(p_nt)
      out <- do.call(c, res)
      attributes(out) <- attributes(res)
    }else{
      out <- p_n
    }
    return(out)
  }
  ## point estimator
  res <- prob_n(data = data)
  res <- as.data.frame(res)
  ## Bootstrap
  prob_boot <- function (data, B, id)
  {
    ids <- unique(data[[id]])
    n <- length(ids)
    res <- matrix(NA, nrow(res), B)
    for (b in 1 : B) {
      bootdata <- NULL
      bids <- sample(ids, replace = TRUE)
      bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
      bootdata <- data[bidxs, ]
      bootres <-  prob_n(bootdata)[, 2]
      res[ , b] <- bootres
    }
    return(res)
  }
  
  ## CI
  if (CI == TRUE) {
    res_boot <- prob_boot(data = data, B = 100, id = "id")
    se <- apply(res_boot, 1, sd)
    
    #############
    ### ? /sqrt(n)
    #############
    res$se <- se  
    ## cloglog transformation
    if (is.null(h)){
      est <- res[[paste0("p", j)]]
    } else {
      est <- res[[paste0("p", h, j)]]
    }
    res$ll <- exp(-exp(log(-log(est)) - qnorm(0.975) * se/
                     (est * log(est))))
    res$ul <- exp(-exp(log(-log(est)) + qnorm(0.975) * se/
                         (est * log(est))))
  }
  return(res)
}
