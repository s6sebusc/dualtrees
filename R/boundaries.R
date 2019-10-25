#' Padded boundary conditions
#' @export
make_square <- function( picture, N, Ny = N , value=min(picture, na.rm=TRUE) ){
    
    nx <- nrow(picture)
    ny <- ncol(picture)
    if( nx > N | ny > Ny ) stop( "give me a bigger square" )
    
    px <- 1:nx + ceiling( (N - nx)/2 ) 
    py <- 1:ny + ceiling( (Ny - ny)/2 ) 
    
    res <- array( dim=c(N,Ny), data=value )
    res[ px, py ] <- picture
    
    l  <- list( res=res, px=px, py=py )
    class( l ) <- "bc_field"
    return( l )
}

#' Reflective boundary conditions
#' @export
put_in_mirror <- function( x, N, Ny=N ){
    x  <- as.matrix(x)
    nx <- nrow( x )
    ny <- ncol( x )
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

#' Periodic boundary conditions
#' @export
period_bc <- function( x, N, Ny=N ){
    x  <- as.matrix( x )
    nx <- nrow( x )
    ny <- ncol( x )
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

