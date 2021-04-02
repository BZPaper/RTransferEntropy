# Function for bootstrapping Renyi transfer entropy under H0 of independence
# between time series x and y. Used internally by transfer_entropy; same
# arguments.
#
bootstrap_renyi <- function(x,
                            lx,
                            y,
                            ly,
                            q,
                            burn = 50) {
  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  # Lead = x
  dteyx <- calc_te_renyi(
    x = bootx,
    lx = lx,
    y = y,
    ly = ly,
    q = q
  )

  # Lead = y
  dtexy <- calc_te_renyi(
    x = booty,
    lx = ly,
    y = x,
    ly = lx,
    q = q
  )

  teboot <- c(dtexy, dteyx)
  names(teboot) <- c("dtexy", "dteyx")

  return(teboot)
}
