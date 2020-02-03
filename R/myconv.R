my_conv <- function( mat, fil, dec=TRUE, odd=FALSE, inverse=FALSE ){

    bc  <- period_bc( mat, N=nrow(mat)+length(fil), Ny=ncol(mat) )

    if( dec | inverse ) mat <- bc$res
    
    if( length( fil ) > nrow( mat ) ) stop( "filter too long" )
    res <- array( dim=dim(mat), data=NA )
    l   <- length(fil)
    lh  <- floor(l/2)
    fil <- c( fil[(lh+1):l], rep( 0, nrow(mat) - l ) , fil[ 1:(lh) ] )
    fil <- matrix( nrow=nrow(mat), ncol=ncol(mat), data=stats::fft( fil ) )
    res <- stats::mvfft( stats::mvfft( mat )*fil, inverse=TRUE ) / nrow(mat)
    
    if( dec | inverse ) res <- res[ bc$px, , drop=FALSE ] 
    if( dec ) res <- decimate( res, odd=odd )
    return( res )
} 
