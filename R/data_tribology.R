#' Tribology data
#'
#' A dataset containing FTIR spectra of engine oils, and
#' HoG (histogram of gradients) feature vectors representing wear scar images
#' after an SRV experiment.
#'
#' @usage
#' data(tribology)
#' @format
#' A list containing the following components:
#' \describe{
#'   \item{FTIR}{a numeric matrix (214 * 1668) containing FTIR spectra}
#'   \item{HOG}{a numeric matrix (214 * 1764) containing HoG features}
#'   \item{duration}{numeric (214) artificial alteration duration of oils}
#'   \item{wavenumbers}{char (1668) column names for FTIR data}
#' }
#' @source The work to obtain this data was funded by the Austrian COMET-Program (project InTribology1, no. 872176) via the
#' Austrian Research Promotion Agency (FFG) and the federal states of Niederösterreich and Vorarlberg and
#' was carried out at the Austrian Excellence Centre of Tribology (AC2T research GmbH).
#' @references <https://www.ac2t.at/en/>
#' @examples
#' # data(tribology)
#'
"tribology"
