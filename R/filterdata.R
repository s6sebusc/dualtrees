#' filterbanks for the dtcwt
#'
#' Some of the filters implemented in the python package dtcwt.
#'
#' @docType data
#'
#'
#' @format A list of high- and low-pass filters for analysis and synthesis
#'
#' @keywords datasets
#'
#' @source dtcwt python package
#' @details The near-sym filterbanks are biorthogonal wavelets used for the first level, they have 13 and 19 taps. The qshift filterbanks, each with 14 taps, are suitable for all higher levels of the dtcwt, as the a- and b-filters for an approximate Hilbert-pair. The naming convention follows the python-package: \tabular{ll}{ h: \tab analysis\cr g:\tab synthesis\cr 0:\tab low-pass\cr 1:\tab high-pass\cr a,b: \tab shifted filters}
#' The \code{b_bp}-versions of the filterbanks contain a second high-pass for the diagonal directions, denoted by 2.
#'
#' @name filterbanks
NULL

#' @rdname filterbanks
"qshift_b"

#' @rdname filterbanks
"qshift_b_bp"

#' @rdname filterbanks
"near_sym_b"

#' @rdname filterbanks
"near_sym_b_bp"
