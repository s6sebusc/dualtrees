#' Various boundary conditions for the 2D wavelet transform.
#'
#' Extend a matrix to the desired size.
#' @param x a real matrix
#' @param N the number of rows of the desired output
#' @param Ny the number of columns of the desired output, defaults to N
#' @param value the value with which the picture is padded by \code{pad}
#' @return a list containing the extended matrix (\code{$res}) and the positions of the original matrix within the extended one (\code{$px} and \code{$py}). 
#' @details \code{pad} pads the fields with a constant value on all sides, be careful what you pick here. \code{put_in_mirror} reflects the input at all edges (with repeated end samples), \code{period_bc} simply repeats the input periodically. In any case, you can retrieve the initial area via \code{bc$res[ bc$px, bc$py ]}. 
#' @note N and Ny must be at least as big as the input.
#' @examples
#' bc <- put_in_mirror( blossom, N=300 )
#' plot( bc )
#' print( range( bc$res[ bc$px, bc$py ] - blossom ) )
#' @name boundaries
NULL

#' @rdname boundaries
#' @export
pad <- function( x, N, Ny = N , value=min(x, na.rm=TRUE) ){
    
    nx <- nrow(x)
    ny <- ncol(x)
    if( nx > N | ny > Ny ) stop( "output dimensions too small" )
    
    px <- 1:nx + ceiling( (N - nx)/2 ) 
    py <- 1:ny + ceiling( (Ny - ny)/2 ) 
    
    res <- array( dim=c(N,Ny), data=value )
    res[ px, py ] <- x
    
    l  <- list( res=res, px=px, py=py )
    class( l ) <- "bc_field"
    return( l )
}

#' @rdname boundaries
#' @export
put_in_mirror <- function( x, N, Ny=N ){
    x  <- as.matrix(x)
    nx <- nrow( x )
    ny <- ncol( x )
    if( nx > N | ny > Ny ) stop( "output dimensions too small" )
    
    if( nx==N & ny==Ny ){
        px <- 1:N
        py <- 1:Ny
    }else{
        dx <- floor( (nx-N)/2 )  
        dy <- floor( (ny-Ny)/2 )
        x1 <- y1 <- 0
        while( ncol( x ) < Ny ){
            y1 <- y1 + ncol(x)
            x <- cbind( x[ ,ncol(x):1, drop=FALSE ], x,  x[ ,ncol(x):1, drop=FALSE ] )
        }
        while( nrow( x ) < N ){
            x1 <- x1 + nrow(x)
            x <- rbind( x[ nrow(x):1, , drop=FALSE], x,  x[ nrow(x):1, , drop=FALSE] )
        }
        x  <- x[ 1:N + x1 + dx , 1:Ny + y1 + dy  , drop=FALSE]
        px <- 1:nx - dx 
        py <- 1:ny - dy 
    }
    l  <- list( res=x, px=px, py=py )
    class( l ) <- "bc_field"
    return( l )

}

#' @rdname boundaries
#' @export
period_bc <- function( x, N, Ny=N ){
    x  <- as.matrix( x )
    nx <- nrow( x )
    ny <- ncol( x )
    if( nx > N | ny > Ny ) stop( "output dimensions too small" )
    
    dx <- floor( (nx-N)/2 )
    dy <- floor( (ny-Ny)/2 )
    x1 <- y1 <- 0

    while( ncol( x ) < Ny ){
        y1 <- y1 + ncol(x)
        x <- cbind( x, x, x )
    }
    while( nrow( x ) < N ){
        x1 <- x1 + nrow(x)
        x <- rbind( x, x, x )
    }
    
    x  <- x[ 1:N + x1 + dx , 1:Ny + y1 + dy, drop=FALSE  ]
    px <- 1:nx - dx 
    py <- 1:ny - dy 
    
    l  <- list( res=x, px=px, py=py )
    class( l ) <- "bc_field"
    return( l )
}


#' smoother borders
#'
#' let a field decrease linearly towards its edges
#' @param x a real matrix
#' @param r a positive integer 
#' @return a matrix of the same dimensions as x
#' @details Values within the field are linearly reduced from their original value to the field minimum, starting \code{r} pixels away from the edge. This enforces truely periodic boundaries and removes sharp edges.
#' @note r must not be larger than \code{min( dim(x) )/2}.
#' @examples
#' image( smooth_borders(blossom, r=64), col=gray.colors(128,0,1) )
#' @export
smooth_borders <- function( x, r ){
    if( r > 0 ){
        mi <- min( x )
        x <- x - mi
        di <- dim( x )
        if( min(di)/2 < r ) stop( "r is too large" )
        mask <- array( dim=di, data=0 )
        mv <- seq( 0,1,,r )
        for( i in 1:r ) mask[ i:(di[1]-i+1), i:(di[2]-i+1) ] <- mv[i]
        x <- x*mask + mi
    }
    return( x )
}

