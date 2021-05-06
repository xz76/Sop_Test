library(matrixStats)
library(truncdist)
library(survival)
library(extraDistr)
library(doMC)

source("sopt_test.R")
source("get_data.R")
source("sop_t.R")

doMC::registerDoMC(cores = parallel::detectCores())

run_sim <- function(N = 1000, n = 20, M0 = 5, M1 = 15, tmat, cid = "cid", 
                    id = "id", group = "group",  j = 2, B = 1000, 
                    method = "linear"){
  out <- foreach(i = seq_len(N)) %dopar%
    ({ 
      tdat1 <- simulate(n = n/2, M0 = M0, M1 = M1)
      tdat1$group <- 0
      tdat2 <- simulate(n = n/2, M0 = M0, M1 = M1)
      tdat2$group <- 1
      tdat <- rbind(tdat1, tdat2)
      res <- try(sopt_test(data = tdat, tmat = tmat, cid = "cid", id = "id",
                           group = "group", j = 2, B = B), silent=TRUE)
      if(class(res)!="try-error"){
        res
      }else{
        NA_real_
      }
    })
  res <- do.call(rbind, out)
  # type1 <- colMeans((res - 0.05) < 0, na.rm = TRUE) 
  # count <- sum(!is.na(res))
  return(res)
}

tmatrix <- trans(state_names = c("health", "illness", "death"),from = c(1, 1, 1, 2, 2),
                 to = c(2, 2, 3, 3, 1))
set.seed(0222)
linear11 <- run_sim(N = 1000, n = 20, M0 = 15, M1 = 30, 
                    tmat = tmatrix, B = 1000)
saveRDS(linear11, "l12.rds")

linear21 <- run_sim(N = 1000, n = 40,  M0 = 15, M1 = 30,  tmat = tmatrix, B = 1000)
saveRDS(linear21, "l22.rds")

linear31 <- run_sim(N = 1000, n = 80, M0 = 15, M1 = 30,  tmat = tmatrix, B = 1000)
saveRDS(linear31, "l32.rds")

# KS12 <- run_sim(N = 10, n = 20, M0 = 5, M1 = 10,  tmat = tmatrix, B = 1000)
# saveRDS(KS11, "k12.rds")
# 
# KS22 <- run_sim(N = 1000, n = 40, M0 = 10, M1 = 30,  tmat = tmatrix, B = 1000)
# saveRDS(KS22, "k22.rds")
# 
# KS32 <- run_sim(N = 1000, n = 80, M0 = 10, M1 = 30,  tmat = tmatrix, B = 1000)
# saveRDS(KS32, "k32.rds")
# 
# KS21 <- run_sim(N = 1000, n = 40, M0 = 5, M1 = 15,  tmat = tmatrix, B = 1000)
# saveRDS(KS21, "k21.rds")
# 
# KS31 <- run_sim(N = 1000, n = 80, M0 = 5, M1 = 15,  tmat = tmatrix, B = 1000)
# saveRDS(KS31, "k31.rds")
# 
# KS12 <- run_sim(N = 1000, n = 20, M0 = 10, M1 = 30,  tmat = tmatrix, B = 1000)
# saveRDS(KS11, "k12.rds")
# 
# KS22 <- run_sim(N = 1000, n = 40, M0 = 10, M1 = 30,  tmat = tmatrix, B = 1000)
# saveRDS(KS22, "k22.rds")
# 
# KS32 <- run_sim(N = 1000, n = 80, M0 = 10, M1 = 30,  tmat = tmatrix, B = 1000)
# saveRDS(KS32, "k32.rds")

# res2 <- readRDS("result/l32.rds")
# mean(res2<0.05, na.rm = TRUE)
# 



# c(sum((tmp$from ==2) & (tmp$to == 2)), sum((tmp$from ==2)))
