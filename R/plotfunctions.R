#' @export
plot.bc_field <- function( bc, lcol="green", col=gray.colors(32,0,1), lty=2, ... ){
    im  <- bc$res
    nx  <- length(bc$px)
    ny  <- length(bc$py)
    x   <- seq(0,1,,nrow(im))
    y   <- seq(0,1,,ncol(im))
    bcx <- x[ bc$px[ c( 1, nx, nx, 1, 1 ) ] ]
    bcy <- y[ bc$py[ c( 1, 1, ny, ny, 1 ) ] ]
    image( im, col=col, ... )
    lines( bcx, bcy, col=lcol, lty=lty )
}

#' plot centre as vectors
#' 
#' display the radial and angular component of the spectrum's centre as arrows.
#'
#' @param uv array of dimension \code{ nx x ny x 2 }, containing the u- and v-component, result of \code{cen2uv}
#' @param z image to show in the background, defaults to the absolute velocity
#' @param col color of the arrows
#' @param zcol color scale for the image
#' @param n number of arrows in one direction
#' @param f factor by which to enlarge the arrows
#' @param code determines the type of arrow, passed to \code{arrows()}
#' @param length length of the arrowhead in inches, does nothing if \code{code=0}
#' @details The pivot of the arrows is at the location to which the u- and v-component belong. By default, no arrowhead is displayed (\code{code=0}) since the egdges detcted by the cdtwt have an orientation but no sign (SW and NE are equivalent). The default size of the arrows is such that a 'velocity' of 1 corresponds to 5% of the shorter image side.
#' @examples
#' uv <- cen2uv( dt2cen( fld2dt( boys ) ) )
#' uvplot( uv, z=boys )
#' @seealso \code{\link{cen2uv}}
#' @export
uvplot <- function(  uv, z=NULL , x=NULL, y=NULL, col="green", zcol=gray.colors(32,0,1), n=42, f=1, code=0, length=.05 ){
    nx <- nrow( uv )
    ny <- ncol( uv )
    if( is.null(z) ) z <- sqrt( uv[,,1]**2 + uv[,,2]**2 )
    if( ( nrow(z) != nx ) | (ncol(z) != ny) ){ 
        stop( "dimensions of uv an z don't match." )
    }
    if( is.null(x) ) x <- seq( 0, 1, , nx )
    if( is.null(y) ) y <- seq( 0, 1, , ny )
    
    f <- f*.05*min( max(x), max(y) )
    
    px  <- seq( 1,nx,,n )
    py  <- seq( 1,ny,,n )
    u  <- c(uv[ px,py,1 ])*f
    v  <- c(uv[ px,py,2 ])*f
    
    xy <- expand.grid(x=x[px], y=y[py])
    image( x,y,z, col=zcol, xaxt="n", yaxt="n", xlab="", ylab="" )
    with( xy, arrows( x0=x-u/2, y0=y-v/2, x1=x+u/2, y1=y+v/2, length=length, col=col, code=code ) )
}
