indval <- function(dat, clas = colnames(dat), spec = rownames(dat), 
                   method = "APCF", rawdata = FALSE, significance = FALSE, 
                   nboot = 999) { 
    
    # Variables description: see the function's previous version
    # The function works by different ways with 0 and NA in matrix:
    #    0 means abcence and involved to calculations
    #    NA means `no data` and skipped of calculations
    
    
    # test & pre-processing unit ----------------------------------------------
    dat <- as.matrix(dat)
    # dat[is.na(dat)] <- 0
    if(!is.numeric(dat)) {show("Your data aren't numeric"); break}
    if(sum(dat, na.rm = TRUE) <= 0) {show("Your data are empty or negative"); break}
    if(length(clas) != ncol(dat)) {show("Classification size error"); break}
    if(length(spec) != nrow(dat)) {show("Species number error"); break}
    if(!(method %in% c("APCF", "ACF", "ASF", "PCF", "PSF"))) {
        show("You selected unknown type of IndVal
try one of this: APCF, ACF, ASF, #PCF, #PSF"); break}
    clas[clas == "NA"] <- NA
    g <- unique(clas) |> sort() |> na.omit()
    
    # Selecting unit ----------------------------------------------------------
    if(method == "APCF") { 
        calcA <- function(d, clas) { 
            res0 <- rep(NA, length(g))
            names(res0) <- g
            res <- rep(NA, length(g))
            names(res) <- g
            for(h in 1:length(g)) {
                res0[h] <- mean(d[clas == g[h]], na.rm = TRUE)
            }
            res0[is.nan(res0)] <- NA
            for(h in g){ 
                res[[h]] <- res0[[h]]/sum(res0, na.rm = TRUE)
            }
            return(res)
        }
        calcB <- function(d, clas) { 
            res <- rep(NA, length(g))
            names(res) <- g
            d[d>0] <- 1
            for(j in g)
                res[[j]] <- mean(d[which(clas == j)], na.rm = TRUE)
            return(res)
        }
    }
    if(method == "ACF") { 
        calcA <- function(d, clas) { 
            res0 <- rep(NA, length(g))
            names(res0) <- g
            res <- rep(NA, length(g))
            names(res) <- g
            for(h in 1:length(g)) {
                res0[h] <- mean(d[clas == g[h]], na.rm = TRUE)
            }
            res0[is.nan(res0)] <- NA
            for(h in g){ 
                res[[h]] <- res0[[h]]/sum(res0, na.rm = TRUE)
            }
            return(res)
        }
        calcB <- function(d, clas) { 
            res <- rep(NA, length(g))
            names(res) <- g
            for(j in g) { 
                e <- d[clas == j] |> na.omit()
                res[[j]] <- 1 - sum(abs(e/sum(e) - 1/length(e)))*0.5
            }
            res[is.nan(res)] <- 0
            return(res)
        }
    }
    if(method == "ASF") { 
        calcA <- function(d, clas) { 
            res0 <- rep(NA, length(g))
            names(res0) <- g
            res <- rep(NA, length(g))
            names(res) <- g
            for(h in g) { 
                res0[[h]] <- mean(d[clas == h], na.rm = TRUE)
            }
            res0[is.nan(res0)] <- NA
            for(h in g) { 
                res[[h]] <- (res0[[h]] - mean(res0[names(res0) != h], 
                                              na.rm = TRUE)) / max(res0, na.rm = TRUE)
            }
            return(res)
        }
        calcB <- function(d, clas) { 
            res <- rep(NA, length(g))
            names(res) <- g
            for(j in g) { 
                e <- d[clas == j] |> na.omit()
                res[[j]] <- 1 - sum(abs(e/sum(e) - 1/length(e)))*0.5
            }
            res[is.nan(res)] <- 0
            return(res)
        }
    }
    
    # computational unit ------------------------------------------------------
    A <- matrix(NA, nrow = length(spec), 
                ncol = length(g))
    colnames(A) <- g
    rownames(A) <- spec
    B <- A
    for(i in 1:nrow(A)) { 
        A[i,] <- calcA(dat[i,], clas)
        B[i,] <- calcB(dat[i,], clas)
    }
    res <- A*B*100
    res[is.nan(res)] <- NA
    S <-  mean(abs(res), na.rm = TRUE)
    
    # permutational test ------------------------------------------------------
    if(significance == TRUE) { 
        library(parallel)
        library(doParallel)
        library(foreach)
        clas2 <- matrix(NA, ncol = length(clas), nrow = nboot)
        for(i in 1:nrow(clas2)) { 
            clas2[i,] <- sample(clas, length(clas))
        }
        myCluster <- makeCluster(detectCores()-1)
        registerDoParallel(myCluster)
        permutated <- foreach(i = 1:nboot, .combine = 'c') %dopar% {
            res2 <- matrix(NA, nrow = length(spec), ncol = length(g))
            for(j in 1:nrow(res2)) {
                res2[j,] <- 100 * 
                    calcA(dat[j,], clas2[i,]) 
                calcB(dat[j,], clas2[i,])
            }
            mean(abs(res2), na.rm = TRUE)
        }
        stopCluster(myCluster)
        S2 <- (S - mean(c(S, permutated), na.rm = TRUE))/sd(c(S, permutated))
    }
    
    # Return results ----------------------------------------------------------
    if(rawdata == FALSE & significance == FALSE) {return(
        list(Method = method, Sharpness = S, indval = res))
    }
    if(rawdata == TRUE  & significance == FALSE) {return(
        list(Method = method, Sharpness = S, indval = res, A = A, B = B))
    }
    if(rawdata == FALSE & significance == TRUE)  {return(
        list(Method = method, Sharpness = S, Significance = S2, indval = res)) 
    }
    if(rawdata == TRUE  & significance == TRUE)  {return(list(Method = method,
                                                              Sharpness = S, Significance = S2, indval = res, A = A, B = B) )
    }
}