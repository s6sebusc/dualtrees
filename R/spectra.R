#' get energy from the dualtree transform
#' 
#' square the wavelet coefficients and apply some form of correction
#' @param pyr array of \code{ J x nx x ny x 6 } complex wavelet coefficients, output of \code{dtcwt(..., dec=FALSE)}
#' @param correct type of correction, either \code{"b"}, \code{"b_bp"} or \code{"fast"}
#' @return an array of the same dimensions as \code{pyr}
#' @details The bias correction matrix should correspond to the filter bank used in the transform. It is computed by brute force for the so-called auto-correlation wavelets, for more details, see Nelson et al 2017 and Eckley et al 2010.
#' @references Eckley, Idris A., Guy P. Nason, and Robert L. Treloar. “Locally Stationary Wavelet Fields with Application to the Modelling and Analysis of Image Texture: Modelling and Analysis of Image Texture.” Journal of the Royal Statistical Society: Series C (Applied Statistics), June 21, 2010, no-no. \url{https://doi.org/10.1111/j.1467-9876.2009.00721.x}.
#' 
#' Nelson, J. D. B., A. J. Gibberd, C. Nafornita, and N. Kingsbury. “The Locally Stationary Dual-Tree Complex Wavelet Model.” Statistics and Computing 28, no. 6 (November 2018): 1139–54. \url{https://doi.org/10.1007/s11222-017-9784-0}.
#'
#' @seealso \code{\link{A}}
#' @export
get_en <- function( pyr, correct="fast", N=ncol( pyr ) ){
    pyr <- Mod( pyr )**2
    if( correct=="fast" ){
        for( j in 1:nrow( pyr ) ) pyr[j,,,] <- pyr[j,,,]/(2**j)
    }
    if( correct=="b_bp" ){
        pyr <- biascor( pyr, A_b_bp[[ paste0( "N", N ) ]] )
    }
    if( correct=="b" ){
        pyr <- biascor( pyr, A_b[[ paste0( "N", N ) ]] )
    }
    return( pyr )
}

#' spectral bias correction, implemented in FORTRAN 
#'
#' Unfold scale and direction into a vector and multiply by A at each grid point
#' @param en \code{ J x nx x ny x 6 } array of squared wavelet coefficients
#' @param a bias correction matrix 
#' @return an array of the same dimensions as en
#' @note this step can introduce negative values into the spctrum which have no intuitive interpretation in terms of energy and need to be dealt with.
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

#' transform a field into an array of spectral energies
#' 
#' Handles the transformation itself, boundary conditions and bias correction and returns the unbiased local wavelet spectrum at each grid-point.
#' @param fld a real matrix
#' @param Nx size to which the field is padded in x-direction
#' @param Ny size to which the field is padded in y-direction
#' @param J number of levels for the decomposition
#' @param mode how to handle the convolutions
#' @param correct how to correct the bias, either "fast", "b" or "b_bp" - any other value results in no correction
#' @param verbose whether or not you want the transform to talk to you
#' @param boundaries how to handle the boundary conditions, either "pad", "mirror" or "periodic"
#' @param fb1 filter bank for level 1
#' @param fb2 filter bank for all further levels
#' @return an array of size \code{J x nx x ny x 6} where \code{dim(fld)=c(nx,ny)}
#' @details The input is blown up to \code{Nx x Ny} and the thrown into \code{dtcwt}. Then the original domain is cut out, the coefficients are squared and the bias is corrected.
#' @examples
#' dt <- fld2dt( boys )
#' par( mfrow=c(2,2), mar=rep(2,4) )
#' for( j in 1:4 ){
#'     image( boys, col=gray.colors(128, 0,1), xaxt="n", yaxt="n" )
#'     for(d in  1:6) contour( dt[j,,,d], levels=quantile(dt[,,,], .995), 
#'                             col=d+1, add=TRUE, lwd=2, drawlabels=FALSE )
#'     title( main=paste0("j=",j) )
#' } 
#' x0  <- seq( .1,.5,,6 )
#' y0  <- rep( 0.01,6 )
#' a   <- .075
#' phi <- seq( 15,,30,6 )*pi/180
#' x1  <- x0 + a*cos( phi )
#' y1  <- y0 + a*sin( phi )
#' rect( min(x0,x1)-.05, min(y0,y1)-.05, 
#'       max(x0,x1)+.05, max(y0,y1), col="black", border=NA )
#' arrows( x0, y0, x1, y1, length=.05, col=2:7, lwd=2, code=3 )

#' @seealso \code{\link{biascor}}
#' @export
fld2dt <- function( fld, Nx=NULL, Ny=NULL, J=NULL, mode=NULL, correct=NULL, verbose=FALSE, boundaries="pad", fb1=near_sym_b_bp, fb2=qshift_b_bp ){
    
    if( is.null( Nx ) ) Nx <- 2**ceiling( log2( max( dim( fld ) )  ) )
    if( is.null( Ny ) ) Ny <- Nx
    if( is.null( J ) ) J <- log2( min(Nx,Ny) ) - 3
    if( is.null(correct) ){
        if( length( fb1 ) > 4   ){
            correct <- "b_bp"
        }else{
            correct <- "b"
        }
    }  
    
    if( boundaries=="mirror" ){
        bc  <- put_in_mirror( fld, N=Nx, Ny=Ny )
    }
    if( boundaries=="pad" ){
        bc  <- pad( fld, N=Nx, Ny=Ny )
    }
    if( boundaries=="periodic" ){
        bc  <- period_bc( fld, N=Nx, Ny=Ny )
    }
    fld <- bc$res
    
    res <- dtcwt( fld, dec=FALSE, mode=mode, J=J, 
                  fb1=fb1, fb2=fb2,
                  verbose=verbose )
    res <- get_en( res[ ,bc$px,bc$py, ], correct=correct, N=2**(J+3) )
    
    return( res )
}

#' spatial mean spectrum
#' 
#' average the output of fld2dt or dtcwt over space
#'
#' @param x either a J x nx x ny x 6 array of energies (output of dtcwt) or a list of complex wavelet coefficients (the output of \code{dtcwt(...,dec=FALSE)})
#' @return a J x 6 matrix of spatially averaged energies
#' @note In the undecimated case, the coefficients are not averaged but summed up and then scaled by the area of the first level. This yields a comparable scale as the undecimated case.
#' @export
dtmean <- function( x ){
    if( is.list(x) ){
        J   <- length(x) - 1
        res <- array( dim=c(J,6) )
        f   <- prod( dim( x[[1]] ) )
        for( j in 1:J ) for( d in 1:6 ) res[ j,d ] <- sum( Mod( x[[j]][,,d] )**2 )/f
    }else{
        res <- apply( x, c(1,4), mean, na.rm=TRUE )
    }
    return( res )
}

 
#' centre of the DT-spectrum
#' 
#' calculate the centre of mass of the local spectra in hexagonal geometry
#' @param pyr a \code{J x nx x ny x 6} array of spectral energies, the output of \code{fld2dt}
#' @return a \code{nx x ny x 3} array where the third dimension denotes degree of anisotropy, angle and central scale, respectively.
#' @details Each of the \code{J x 6} spectral values is assigned a coordinate in 3D space with \code{x(d,j)=cos(60*(d-1))}, \code{y(d,j)=sin(60*(d-1))},\code{z(d,j)=j+1}. Then the centre of mass in this space is calculated, the spectral values being the masses at each vertex. The x- and y-cooridnate are then transform into a radius \code{rho=sqrt(x^2+y^2)} and and angle \code{phi=15+0.5*atan2(y,x)}. \code{rho} measures the degree of anisotropy at each pixel, \code{phi} the orientation of edges in the image, and the third coordinate, \code{z}, the central scale.
#' @note Since the centre of mass is not defined for negative mass, any values below zero are removed at this point. 
#' @examples
#' dt <- fld2dt(blossom)
#' ce <- dt2cen(dt)
#' image( ce[,,3], col=gray.colors(32, 0, 1) )
#' @useDynLib dualtrees
#' @export
dt2cen <- function( pyr ){
    if( length( dim(pyr) ) == 2 ) pyr <- array( dim=c( nrow(pyr),1,1,ncol(pyr) ), data=pyr )
    pyr[pyr<0] <- 0
    nx  <- as.integer( dim(pyr)[2] )
    ny  <- as.integer( dim(pyr)[3] )
    nl  <- as.integer( dim(pyr)[1] )
    tmp <- array( dim=c( nx, ny, 3 ), data=0 ) 
    tmp <- .Fortran( "getcen", pyramid=pyr, nx=nx, ny=ny, nl=nl, res=tmp, PACKAGE="dualtrees" )$res
    res <- tmp
    res[,,1]   <- sqrt( tmp[,,1]**2 + tmp[,,2]**2 )
    phi        <- ( atan2( tmp[,,2],tmp[,,1] )*180/pi )/2 + 15
    phi[phi<0] <- 180 + phi[phi<0]
    res[,,2]   <- phi
    return( res )
}

#' centre to vector
#'
#'transforms the angle and radius component of the spectral centre into vector components for visualization
#' @param cen an \code{nx x ny x 3} array, the output of \code{dt2cen}
#' @return an \code{nx x ny x 2} array, the two matrices representing the u- and v-component of the vector field. 
#' @details The resulting vector field represents the degree of anisotropy and dominant direction, averaged over all scales. In principle, large and small edges might compensate one another, but the result typically looks as one would intuitively expect.
#' @examples
#' dt <- fld2dt(blossom)
#' ce <- dt2cen(dt)
#' uv <- cen2uv(ce)
#' uvplot( uv, z=blossom )
#' @seealso \code{\link{dt2cen}}, \code{\link{uvplot}}
#' @export
cen2uv <- function( cen ){
    if( is.null( dim( cen ) ) ) cen <- array( dim=c(1,1,3), data=cen )
    uv      <- array( dim=c( nrow(cen), ncol(cen), 2 ), data=NA )
    uv[,,1] <- cen[,,1]*cos( cen[,,2]*pi/180 )
    uv[,,2] <- cen[,,1]*sin( cen[,,2]*pi/180 )
    return( uv )
}
