#' Column-convolutions
#'
#' This function convolves the columns of a matrix mat with a filter fil.
#' @param mat a matrix
#' @param fil the filter to convolve the columns with
#' @param dec if \code{TRUE}, every second row is discarded after the convolution
#' @param mode how to actually do the convolutions, must be either \code{"direct"} or \code{"FFT"}
#' @param odd if \code{TRUE}, the first row is discarded, otherwise the second row is.
#' @param how to handle the boundaries, does nothing if \code{mode="FFT"}
#' @keywords convolution, wavelets
#' @details This functions does all of the actual computations inside the wavelet transform. The direct mode uses \code{filter(...)} and can handle any field size you like. It is supposedly faster when the filters are short, i.e., in the decimated case. The FFT-version really only works when the input dimensions are whole powers of two and the filter is not longer than the columns of the matrix.
#' @examples
#' dboysdy <- my_conv( boys, c(-1,1), dec=FALSE )
#' dboysdx <- t( my_conv( t(boys), c(-1,1), dec=FALSE ) )
#' par( mfrow=c(1,2) )
#' image( dboysdx, col=gray.colors(32) )
#' image( dboysdy, col=gray.colors(32) )
#' @export
my_conv <- function( mat, fil, dec=TRUE, mode="direct", odd=FALSE, boundaries="periodic" ){
    nx  <- nrow(mat)
    if( !( mode %in% c( "FFT", "direct" ) ) ){ 
        stop( paste0( "unknown mode '", mode, "', chose 'FFT' or 'direct', friend." ) )
    }
    
    if( mode == "FFT" ){
        if( length( fil ) > nrow( mat ) ) stop( "filter too long, try mode = 'direct'" )
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
