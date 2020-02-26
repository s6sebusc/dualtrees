#' Bias correction matrices
#'
#' Matrices needed in order to eliminate the effect of "spectral leakage" for the local dtcwt-spectra.
#'
#' @docType data
#'
#'
#' @format A list with entries N1024, N512, N256, N128, N64, N32, each containing the bias correction matrix of appropriate size.
#'
#' @keywords datasets
#'
#' @details As descibed in Nelson et al. (2018), the squared coefficients of the undecimated dtcwt (the local wavelet spectrum) can be used to obtain an estimate of local correlations in space. This estimator has a bias which mostly consists of an over-emphasis on the largest scales. It can be removed by multiplying each local spectrum with a matrix which depends only on the choice of wavelet and the dimensions of the field.
#' @source Calculated by brute force.
#' @references Nelson, J. D. B., A. J. Gibberd, C. Nafornita, and N. Kingsbury (2018) <doi:10.1007/s11222-017-9784-0>.
#'
#' @examples
#' image( A_b_bp$N512 )
#' @name A
NULL


#' @rdname A
"A_b_bp"  

#' @rdname A
"A_b"  
