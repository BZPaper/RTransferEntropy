#' Function that generates clusters and calculates frequencies for each cluster. 
#' Time series x and y are assumed to follow a Markov process of order k and j, 
#' respectively.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#' @param prog if TRUE, x(k+1) 
#' @param lx x(k)
#' @param ly y(j)
#'
#' @return returns a list with clusters and associated frequencies
#' @export
#'
#' @examples
#' 
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