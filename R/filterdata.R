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
#' @source 'dtcwt' python package (\url{https://github.com/rjw57/dtcwt})
#' @details The near-sym filterbanks are biorthogonal wavelets used for the first level, they have 13 and 19 taps. The qshift filterbanks, each with 14 taps, are suitable for all higher levels of the dtcwt. The a- and b-filters form an approximate Hilbert-pair. The naming convention follows the python-package: \tabular{ll}{ h: \tab analysis\cr g:\tab synthesis\cr 0:\tab low-pass\cr 1:\tab high-pass\cr a,b: \tab shifted filters}
#' The \code{b_bp}-versions of the filterbanks contain a second high-pass for the diagonal directions, denoted by 2. They allow for better directional selectivity but prohibit perfect reconstruction.
#' @references Selesnick, I.W., R.G. Baraniuk, and N.C. Kingsbury (2005) <doi:10.1109/MSP.2005.1550194>
#'
#' Kingsbury, N. (2006) <https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7071567&isnumber=7065146>
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
