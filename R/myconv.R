#' Column-convolutions
#'
#' This function convolves the columns of a matrix mat with a filter fil.
#' @param mat a matrix
#' @param fil the filter to convolve the columns with
#' @param dec if \code{TRUE}, every second row is discarded after the convolution
#' @param odd if \code{TRUE}, the first row is discarded, otherwise the second row is.
#' @param inverse set \code{TRUE} for the inverse transformation
#' @keywords convolution, wavelets
#' @details This functions does all of the actual computations inside the wavelet transform via FFT. If dec=TRUE or inverse=TRUE, the filters may be longer than the data. In these cases, the data is periodically extended to exceed the filter size.
#' @return a matrix with as many columns and either the same (\code{dec=FALSE}) or half (\code{dec=TRUE}) the number of rows as mat 
#' @examples
#' dboysdy <- my_conv( blossom, c(-1,1), dec=FALSE )
#' dboysdx <- t( my_conv( t(blossom), c(-1,1), dec=FALSE ) )
#' par( mfrow=c(1,2) )
#' image( dboysdx, col=gray.colors(32) )
#' image( dboysdy, col=gray.colors(32) )
#' @export
my_conv <- function( mat, fil, dec=TRUE, odd=FALSE, inverse=FALSE ){

    bc  <- period_bc( mat, N=nrow(mat)+length(fil), Ny=ncol(mat) )

    if( dec | inverse ) mat <- bc$res
    
    if( length( fil ) > nrow( mat ) ) stop( "filter too long" )
    res <- array( dim=dim(mat), data=NA )
    l   <- length(fil)
    lh  <- floor(l/2)
    fil <- c( fil[(lh+1):l], rep( 0, nrow(mat) - l ) , fil[ 1:(lh) ] )
    fil <- matrix( nrow=nrow(mat), ncol=ncol(mat), data=fft( fil ) )
    res <- mvfft( mvfft( mat )*fil, inverse=TRUE ) / nrow(mat)
    
    if( dec | inverse ) res <- res[ bc$px, , drop=FALSE ] 
    if( dec ) res <- decimate( res, odd=odd )
    return( res )
} 
