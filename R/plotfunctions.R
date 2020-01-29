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

#' @export
plot.dtana <- function( x ){   
    J <- nrow(x$ms)
    with( x, {
    par( mfrow=c(2,3), mar=c(3,3,3,1) )
    plot( rmi, rhist, type="h", main="histogram of radii" )
    plot( pmi, phist, type="h", main="histogram of angles" )
    plot( zmi, zhist, type="h", main="histogram of central scales" )
    image( 1:J, 1:6, ms, yaxt="n", main="mean spectrum (white = large)" )
    axis( 2, at=1:6, labels=paste0(seq(15,,30,6),"°") )
    })
    if( !is.null(x$uv) ){ 
        uvplot( x )
        title( main="direction and anisotropy" )
        image( x$cen[,,3], xaxt="n", yaxt="n", zlim=c(1,J), main="central scales (white = large)" )
    }
}

#' plot centre as vectors
#' 
#' display the radial and angular component of the spectrum's centre as arrows.
#'
#' @param uv an array of dimension \code{ nx x ny x 2 }, containing the u- and v-component, result of \code{cen2uv}
#' @param z image to show in the background, defaults to \code{sqrt(x^2+y^2)}
#' @param x,y optional x- and y-coordinates for the plot, must match the dimensions of \code{z}
#' @param col color of the arrows
#' @param zcol color scale for the image
#' @param n number of arrows in one direction
#' @param f factor by which to enlarge the arrows
#' @param length length of the arrowhead in inches
#' @param ... further arguments passed to \code{image}
#' @details The pivot of the arrows is at the location to which the u- and v-component belong. No arrowhead is displayed since the egdges detcted by the cdtwt have an orientation but no sign. The default size of the arrows is such that a 'velocity' of 1 corresponds to 5\% of the shorter image side.
#' @examples
#' uv <- cen2uv( dt2cen( fld2dt( blossom ) ) )
#' uvplot( uv, z=blossom )
#' @seealso \code{\link{cen2uv}}
#' @export
uvplot <- function(  uv, z=NULL , x=NULL, y=NULL, col="green", zcol=gray.colors(32,0,1), n=42, f=1, length=.05, ... ){
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
    image( x,y,z, col=zcol, xaxt="n", yaxt="n", xlab="", ylab="", ... )
    with( xy, arrows( x0=x-u/2, y0=y-v/2, x1=x+u/2, y1=y+v/2, length=length, col=col, code=0 ) )
}
