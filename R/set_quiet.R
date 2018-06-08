#' Set the quiet-parameter for all RTransferEntropy Calls
#'
#' @param quiet if FALSE, the functions will give feedback on the progress
#'
#' @return nothing
#' @export
#'
#' @examples
#' # see ?transfer_entropy
set_quiet <- function(quiet) {
  options("RTransferEntropy::quiet" = quiet)
  invisible(NULL)
}
