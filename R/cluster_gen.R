#' Function that generates clusters of states and calculates frequencies for
#' each cluster.
#'
#' @param x a vector of coded values
#' @param y a vector of coded values
#' @param prog if TRUE, lag of x (Markov order) is increased by one
#' @param lx Markov order of x
#' @param ly Markov order of y
#'
#' @return returns a list with clusters and associated frequencies
#' @keywords internal
#' @export
#'
#' @examples
#'
cluster_gen <- function(x,
                        lx,
                        y = NULL,
                        ly = NULL,
                        prog = TRUE) {

  nclust <- length(x) - max(lx, ly)
  clustlist <- list()

  x_pre  <- if (is.null(x)) 0 else max(ly - lx, 0)
  x_post <- lx + if (is.null(x)) 0 else max(ly - lx, 0)
  x_post <- x_post - if (!prog) 1 else 0
  y_pre  <- max(lx - ly, 0)
  y_post <- max(lx, ly) - 1

  x_selector <- x_pre:x_post
  y_selector <- y_pre:y_post

  if (is.null(y)) {
    for (i in 1:nclust) {
      clustlist[[i]] <- x[x_selector + i]
    }
  } else {
    for (i in 1:nclust) {
      clustlist[[i]] <- c(x[x_selector + i],
                          y[y_selector + i])
    }
  }

  numclust <- unlist(lapply(clustlist, function(x) paste(x, collapse = " ")))
  freq <- table(numclust)/length(numclust)

  return(list(cluster = numclust,
              frequency = freq))
}
