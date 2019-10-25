#' get energy from the dualtree transform
#' @export
get_en <- function( pyr, correct="fast" ){
    pyr <- Mod( pyr )**2
    if( correct=="fast" ){
        for( j in 1:nrow( pyr ) ) pyr[j,,,] <- pyr[j,,,]/(2**j)
    }
    if( correct=="b_bp" ){
        N <- 2**( nrow( pyr ) + 3 )
        infile <- paste0( "/user/s6sebusc/wavelets_verification/new_wavelets/Amats/near_sym_b_bp", N, ".rdata" )
        load( infile )
        pyr <- biascor( pyr, A )
    }
    return( pyr )
}

#' spectral bias correction, implemented in FORTRAN 
#' @useDynLib dualtrees
#' @export
biascor <- function( en, a ){
    nx <- as.integer( dim(en)[2] )
    ny <- as.integer( dim(en)[3] )
    nl <- as.integer( dim(en)[1] )
    nd <- as.integer( dim(en)[4] )
    res <- .Fortran( "biascor", sp=en, a=a, nx=nx, ny=ny, nl=nl, nd=nd, res=en, PACKAGE="dualtrees" )$res
    return( res )
}

#' transform a field, handle boundary conditions, return energy
#' @export
fld2dt <- function( fld, Nx=NULL, Ny=NULL, J=NULL, mode=NULL, correct="b_bp", verbose=FALSE, boundary="pad" ){
    
    if( is.null( Nx ) ) Nx <- 2**ceiling( log2( max( dim( fld ) )  ) )
    if( is.null( Ny ) ) Ny <- Nx
    if( is.null( J ) ) J <- log2( min(Nx,Ny) ) - 3
    
    if( boundary=="mirror" ){
        bc  <- put_in_mirror( fld, N=Nx, Ny=Ny )
    }
    if( boundary=="pad" ){
        bc  <- make_square( fld, N=Nx, Ny=Ny )
    }
    fld <- bc$res
    
    res <- dtcwt( fld, dec=FALSE, mode=mode, J=J, 
                  fb1=near_sym_b_bp, fb2=qshift_b_bp,
                  verbose=verbose )
    res <- get_en( res, correct=correct )[ ,bc$px,bc$py, ]
    
    return( res )
}
 
#' get the centre of the DT-spectrum
#' @useDynLib dualtrees
#' @export
dt2cen <- function( pyr ){
    if( length( dim(pyr) ) == 2 ) pyr <- array( dim=c( nrow(pyr),1,1,ncol(pyr) ), data=pyr )
    pyr <- aperm( pyr, c(2,3,1,4) )
    pyr[pyr<0] <- 0
    nx  <- as.integer( nrow( pyr ) )
    ny  <- as.integer( ncol( pyr ) )
    nl  <- as.integer( dim(pyr)[3] )
    res <- array( dim=c( nx, ny, 3 ), data=0 ) 
    res <- .Fortran( "getcen", pyramid=pyr, nx=nx, ny=ny, nl=nl, res=res, PACKAGE="dualtrees" )$res
    return( res )
}

#' transform angle and anisotropy into vector components for plotting
#' @export
cen2uv <- function( cen ){
    uv      <- array( dim=c( nrow(cen), ncol(cen), 2 ), data=NA )
    uv[,,1] <- cen[,,1]*cos( cen[,,2]*pi/180 )
    uv[,,2] <- cen[,,1]*sin( cen[,,2]*pi/180 )
    return( uv )
}
