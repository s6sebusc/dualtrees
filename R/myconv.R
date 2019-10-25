#' Convolve the columns of a matrix in a varitey of ways
#'
#' This function convolves the columns of a matrix mat with a filter fil.
#' @param mat a matrix
#' @param fil the filter to convolve the columns with
#' @keywords convolution, wavelets
#' @export
#' @examples
#' require( fields )
#' data( lennon )
#' my_conv( lennon, c(-1,1) )

my_conv <- function( mat, fil, dec=TRUE, mode="direct", odd=FALSE, boundaries="periodic" ){
    nx  <- nrow(mat)
    
    if( !( mode %in% c( "FFT", "direct" ) ) ){ 
        stop( paste0( "unknown mode '", mode, "', chose 'FFT' or 'direct', friend." ) )
    }
    
    if( mode == "FFT" ){
        res <- array( dim=dim(mat), data=NA )
        l   <- length(fil)
        lh  <- floor(l/2)
        fil <- c( fil[(lh+1):l], rep( 0, nrow(mat) - l ) , fil[ 1:(lh) ] )
        fil <- fft( fil )
        fil <- matrix( nrow=nrow( mat ), ncol=ncol( mat ), data=fil )
        res <- mvfft( mvfft( mat )*fil, inverse=TRUE )
        res <- res / nrow( mat )
    }
    if( mode == "direct" ){
        if( !( boundaries %in% c( "periodic", "reflective" ) ) ){ 
            stop( paste0( "unknown boundaries '", boundaries, "', chose 'periodic' or 'reflective', friend." ) )
        }
        if( boundaries == "periodic" ){
            bc  <- period_bc( mat, N=nrow(mat)+length(fil), Ny=ncol(mat) )
        }
        if( boundaries == "reflective" ){
            bc  <- put_in_mirror( mat, N=nrow(mat)+length(fil), Ny=ncol(mat) )
        }
        mat <- bc$res
        res <- filter( mat, fil, method="convolution", circular=TRUE )
        res <- res[ bc$px, , drop=FALSE ] 
    }
    
    if( dec ) res <- decimate( res, odd=odd )
    return( res )
} 
