# Function for bootstrapping Shannon transfer entropy under H0 of independence
# between time series x and y. Used internally by transfer_entropy; same
# arguments.
#
bootstrap_shannon <- function(x,
                              lx,
                              y,
                              ly,
                              burn = 50) {

  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  # Lead = x
  dteyx <- calc_te_shannon(x = bootx,
                           lx = lx,
                           y = y,
                           ly = ly)

  # Lead = y
  dtexy <- calc_te_shannon(x = booty,
                           lx = ly,
                           y = x,
                           ly = lx)

  teboot <- c(dteyx, dtexy)
  names(teboot) <- c("dteyx", "dtexy")

  return(teboot)
}
