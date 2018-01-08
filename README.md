# Transfer-Entropy
Code to implement transfer entropy (Shannon)

# Usage
``` r
library(RTransferEntropy)
set.seed(42)

dt <- data.table(series = rnorm(10000))
dt[, sample := code_sample(series)]
dt
#>            series sample
#>     1:  1.3709584      2
#>     2: -0.5646982      2
#>     3:  0.3631284      2
#>     4:  0.6328626      2
#>     5:  0.4042683      2
#>    ---                  
#>  9996: -0.5653259      2
#>  9997: -1.4819016      2
#>  9998: -0.2989931      2
#>  9999:  0.2146498      2
#> 10000:  0.7732340      2
```
# Example using simulated data 
Simulate a simple model to obtain two time series that are not independent (see simulation study in Dimpfl and Peter (2013)),
i.e. one time series is lag of the other plus noise. In this case, one expects significant information flow from x to y 
and none from y to x.

``` r
set.seed(20170108)
n <- 100000
x <- rep(0, n + 1)
y <- rep(0, n + 1)

for (i in seq(n)) {
  x[i + 1] <- 0.2 * x[i] + rnorm(1, 0, 2)
  y[i + 1] <- x[i] + rnorm(1, 0, 2)
}

x <- x[-1]
y <- y[-1]
```
# For David
``` r
code_sample <- function(x,
                        type = "quantiles",
                        quantiles = c(5, 95),
                        bins = NULL,
                        limits = NULL,
                        scale = 1e10) {

  if (type %in% c("bins", "limits")) {
    UB <- max(x)
    LB <- min(x)

    # find the respective OSeq for the time series
    if (type == "bins") {
      OSeq <- LB + ((UB - LB) / bins) * (0:(bins))
    } else {
      limits <- sort(limits)
      OSeq <- c(LB, limits, UB)
    }

    OSeq[length(OSeq)] <- UB + 1

    for (j in 1:(length(OSeq) - 1)) {
      x[x >= OSeq[j] & x < OSeq[j + i]] <- j * scale
    }
  } else if (type == "quantiles") {
    Qtl <- quantile(x, type = 8, probs = quantiles/100)
    Qtl <- c(Qtl, max(x))
    Qlength <- length(Qtl)

    for (j in 1:Qlength){
      x[x <= Qtl[j]] <- j * scale
    }
  }
  x <- x / scale

  return(x)
}

cluster_gen <- function(x,
                        y = NULL,
                        lx,
                        ly =NULL,
                        prog = TRUE) {

  n <- length(x)

  if (!is.null(ly) && ly>lx) {
    dclust <- ly + 1
    difflag <- ly - lx
  } else {
    dclust <- lx + 1
  }

  nclust <- (n - dclust) + 1
  clustlist <- list()

  if (is.null(y)) {
    if (prog) {
      for (i in 1:nclust) {
        clustlist[[i]] <- x[i:(i + lx)]
      }
    } else {
      for (i in 1:nclust) {
        clustlist[[i]] <- x[i:(i + lx - 1)]
      }
    }
  } else {
    if (lx >= ly) {
      if (prog) {
        for (i in 1:nclust) {
          clustlist[[i]] <- c(x[i:(i + lx)], y[(i + lx - ly):(i + lx - 1)])
        }
      } else {
        for (i in 1:nclust) {
          clustlist[[i]] <- c(x[i:(i + lx - 1)], y[(i + lx - ly):(i + lx - 1)])
        }
      }
    } else {
      if (prog) {
        for (i in 1:nclust) {
          clustlist[[i]] <- c(x[(i + difflag):(i + difflag + lx)],
                              y[i:(i + ly - 1)])
        }
      } else {
        for (i in 1:nclust) {
          clustlist[[i]] <- c(x[(i + difflag):(i + difflag + lx - 1)],
                              y[i:(i + ly - 1)])
        }
      }
    }
  }

  numclust <- unlist(lapply(clustlist, function(x) paste(x, collapse = "")))
  freq <- table(numclust)/n

  return(list(cluster = numclust,
              frequency = freq))
}

transfer_entropy <- function(x, y, lx, ly){

  # Frequencies
  #------------------------------
  # x(k+1) and y(j)
  k1_j <- cluster_gen(x, y, lx = lx, ly = ly)$frequency
  nck1_j <- length(k1_j)

  # x(k+1)
  k1 <- cluster_gen(x, lx = lx)$frequency
  nck1 <- length(k1)

  # x(k) and y(j)
  k_j <- cluster_gen(x, y, lx = lx, ly = ly, prog = FALSE)$frequency
  nck_j <- length(k_j)

  # x(k)
  k <- cluster_gen(x, lx = lx, prog = FALSE)$frequency
  nck <- length(k)

  # Transfer entropy
  entropy <- numeric(nck1_j)
  for(i in 1:nck1_j){
    p1 <- k1[substr(names(k1_j[i]), 1, (lx + 1))]
    p2 <- k_j[paste(unlist(strsplit(names(k1_j[i]), split = NULL))[-(lx + 1)],
                    collapse = "")]
    p3 <- k[substr(names(k1_j[i]), 1, lx)]
    entropy[i] <- k1_j[i] * log2((k1_j[i] * p3) / (p2 * p1))
  }

  return(list(transentropy = sum(entropy),
              numclassk1_j = nck1_j,
              numclassk1 = nck1,
              numclassk_j = nck_j,
              numclassk = nck))
}

shuffled_transfer_entropy <- function(nreps = 2,
                                      shuffles = 6,
                                      diff = TRUE,
                                      x,
                                      lx,
                                      y,
                                      ly,
                                      ncores = parallel::detectCores() - 1) {

  n <- length(x)

  cl <- parallel::makeCluster(ncores)
  on.exit({
    parallel::stopCluster(cl)
  })

  parallel::clusterExport(cl, "nreps")

  shuffle <- parallel::parLapply(cl, seq(shuffles), function(i) {
    res <- replicate(nreps,
                     transfer_entropy(x = x,
                                      y = sample(y, n, replace = TRUE),
                                      lx = lx, ly = ly)$transentropy)
    return(res)
  })

  ste <- mean(as.numeric(as.character(unlist(shuffle))))

  if (diff) {
    te <- transfer_entropy(x = x, y = y, lx = lx, ly = ly)$transentropy - ste
  } else {
    te <- ste
  }

  return(te)
}

x <- code_sample(x)
y <- code_sample(y)
shuffled_transfer_entropy(x, lx = 1, y, ly = 1)
```
