#' Bias correction matrices
#'
#' Matrices supposedly needed in order to eliminate the effect of "spectral leakage" for the local dtcwt-spectra. Used by biascor( ... ).
#'
#' @docType data
#'
#'
#' @format A list with entries N512, N256, N128, N64, N32, each containing the bias correction matrix of appropriate size
#'
#' @keywords datasets
#'
#' @source Calculated by hand via /user/s6sebusc/wavelets_verification/general_scripts/Amats_cdtwt.r
#'
#' @examples
#' image( A_b_bp$N512 )
#' @name A
NULL


#' @rdname A
"A_b_bp"  

#' @rdname A
"A_b"  
