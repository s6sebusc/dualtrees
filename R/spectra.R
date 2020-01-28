#' get energy from the dualtree transform
#' 
#' square the wavelet coefficients and apply some form of correction
#' @param pyr array of \code{ J x nx x ny x 6 } complex wavelet coefficients, output of \code{dtcwt(..., dec=FALSE)}
#' @param correct type of correction, either \code{"b"} or \code{"b_bp"}, any other value results in no correction at all.
#' @return an array of the same dimensions as \code{pyr}
#' @details The bias correction matrix should correspond to the filter bank used in the transform. It is computed by brute force for the so-called auto-correlation wavelets, for more details, see Nelson et al 2017 and Eckley et al 2010.
#' @references Eckley, Idris A., Guy P. Nason, and Robert L. Treloar. “Locally Stationary Wavelet Fields with Application to the Modelling and Analysis of Image Texture: Modelling and Analysis of Image Texture.” Journal of the Royal Statistical Society: Series C (Applied Statistics), June 21, 2010, no-no. \url{https://doi.org/10.1111/j.1467-9876.2009.00721.x}.
#' 
#' Nelson, J. D. B., A. J. Gibberd, C. Nafornita, and N. Kingsbury. “The Locally Stationary Dual-Tree Complex Wavelet Model.” Statistics and Computing 28, no. 6 (November 2018): 1139–54. \url{https://doi.org/10.1007/s11222-017-9784-0}.
#'
#' @seealso \code{\link{A}}
#' @export
get_en <- function( pyr, correct="none", N=ncol( pyr ) ){
    pyr <- Mod( pyr )**2
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
    if( is.null(a) ) stop( "The bias correction matrix you need does not exist." )
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
#' @param correct how to correct the bias, either "b" or "b_bp" - any other value results in no correction
#' @param rsm number of pixels to be linearly smoothed along each edge before applying the boundary conditions (see \code{\link{smooth_borders}}).
#' @param verbose whether or not you want the transform to talk to you
#' @param boundaries how to handle the boundary conditions, either "pad", "mirror" or "periodic"
#' @param fb1 filter bank for level 1
#' @param fb2 filter bank for all further levels
#' @return an array of size \code{J x nx x ny x 6} where \code{dim(fld)=c(nx,ny)}
#' @details The input is blown up to \code{Nx x Ny} and the thrown into \code{dtcwt}. Then the original domain is cut out, the coefficients are squared and the bias is corrected.
#' @examples
#' dt <- fld2dt( blossom )
#' par( mfrow=c(2,2), mar=rep(2,4) )
#' for( j in 1:4 ){
#'     image( blossom, col=gray.colors(128, 0,1), xaxt="n", yaxt="n" )
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
fld2dt <- function( fld, Nx=NULL, Ny=NULL, J=NULL, correct=NULL, rsm=0, verbose=FALSE, boundaries="pad", fb1=near_sym_b_bp, fb2=qshift_b_bp ){
    
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
    fld <- smooth_borders( fld, rsm )
    
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
    
    res <- dtcwt( fld, dec=FALSE, J=J, 
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
#' @param mask a \code{ nx x ny } array of logical values
#' @return a \code{nx x ny x 3} array where the third dimension denotes degree of anisotropy, angle and central scale, respectively.
#' @details Each of the \code{J x 6} spectral values is assigned a coordinate in 3D space with \code{x(d,j)=cos(60*(d-1))}, \code{y(d,j)=sin(60*(d-1))},\code{z(d,j)=j+1}. Then the centre of mass in this space is calculated, the spectral values being the masses at each vertex. The x- and y-cooridnate are then transform into a radius \code{rho=sqrt(x^2+y^2)} and and angle \code{phi=15+0.5*atan2(y,x)}. \code{rho} measures the degree of anisotropy at each pixel, \code{phi} the orientation of edges in the image, and the third coordinate, \code{z}, the central scale. If a \code{mask} is provided, values where \code{mask==TRUE} are set to \code{NA}.
#' @note Since the centre of mass is not defined for negative mass, any values below zero are removed at this point. 
#' @examples
#' dt <- fld2dt(blossom)
#' ce <- dt2cen(dt)
#' image( ce[,,3], col=gray.colors(32, 0, 1) )
#' @useDynLib dualtrees
#' @export
dt2cen <- function( pyr, mask=NULL ){
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
    if( !is.null(mask) ) for( i in 1:3 ) res[ ,,i ][mask] <- NA 
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

#' xy <-> cen
#'
#' Translate the centre of mass back and forth between polar and cartesian coordinates
#' @param cen the centre of mass of a wavelet spectrum (rho, phi, z), output of \code{dt2cen}
#' @param xy the centre of mass in cartesian coordinates (x, y, z), output of cen2xy
#' @details \code{dt2cen} represents the sepctrum's centre in cylinder coordinates because that is more intuitive than the x-y-z position within the hexagonal geometry. If you want to compare two spectra, it makes more sense to consider their distance in terms of x1-x2, y1-y2 since the difference in angle is only meaningful for reasonably large radii. These functions allow you to translate back and forth between the two coordinate systems. 
#' @note \code{cen2xy} is not the same thing as \code{cen2uv} !
#' @seealso \code{\link{dt2cen}}, \code{\link{cen2uv}}
#' @name cen_xy
NULL

#' @rdname cen_xy
#' @export
cen2xy <- function( cen ){
    if( is.null( dim( cen ) ) ) cen <- array( dim=c(1,1,length(cen)), data=cen )
    rho <- cen[,,1,drop=FALSE]
    phi <- cen[,,2,drop=FALSE]
    phi[which(phi>90)] <- phi[which(phi>90)] - 180
    phi <- 2*(phi - 15 )
    phi <- phi*pi/180
    cen[ ,,1 ] <- rho*cos( phi )
    cen[ ,,2 ] <- rho*sin( phi )
    return( cen )
}

#' @rdname cen_xy
#' @export
xy2cen <- function( xy ){
    if( is.null( dim( xy ) ) ) xy <- array( dim=c(1,1,length(xy)), data=xy )

    rho <- sqrt( xy[,,1,drop=FALSE]**2 + xy[,,2,drop=FALSE]**2 )
    phi <- (atan2(xy[, , 2,drop=FALSE], xy[, , 1,drop=FALSE]) * 180/pi)/2 + 15
    phi[which(phi < 0)] <- 180 + phi[which(phi < 0)]
    xy[ ,,1 ] <- rho
    xy[ ,,2 ] <- phi
    return(xy)
}

#' complete dual-tree analysis
#'
#' Transform an input field and calculate a number of summary quantities like the histograms for the three central components, the mean spectrum and the centre of the mean spectrum.
#' @param fld a real matrix to be analysed
#' @param neg_rm whether or not negative values are immediately removed from the local spectra - recommended.
#' @param rbr,pbr,zbr either the number of breaks to use for the histograms of rho, phi and z or vectors containing those breaks. 
#' @param mask a boolean matrix telling the function which pixels to remove from the analyis (see \code{\link{dt2cen}})
#' @param return_fields whether or not you want the function to return the 2D fields of rho, phi, z, u, v and fld itself. 
#' @param ... further parameters passed to \code{\link{fld2dt}}, including wavelet selection and boundary conditions
#' @return and object of class \code{"dtana"}, containing the following: 
#' \tabular{ll}{ \code{ms}    \tab the mean spectrum \cr
#'               \code{ce}    \tab the centre of mass of \code{ms} (rho, phi, z)\cr
#'               \code{rhist} \tab the histogram of radii\cr
#'               \code{phist} \tab the histogram of angles\cr
#'               \code{zhist} \tab the histogram of central scales\cr
#'               \code{rmi}   \tab the bin centres for \code{rhist}\cr
#'               \code{pmi}   \tab the bin centres for \code{phist}\cr
#'               \code{zmi}   \tab the bin centres for \code{zhist}
#'             }
#' if you set \code{ return_fields=TRUE } (the default), the resulf also contains
#' \tabular{ll}{ \code{cen}    \tab the output of \code{dt2cen} \cr
#'               \code{dt}     \tab the output of \code{fld2dt}\cr
#'               \code{uv}     \tab the centre converted into cathesian coordinates\cr
#'               \code{fld}    \tab the input field
#'             }
#' @details The input is transformed via \code{fld2dt( fld, ... )} and the plugged into \code{dt2cen}. The results are summarized as histograms of radii, angles and central scales. By default, 33 equally spaced breaks between the theoretical extrema of these quantities are used. The default breaks for rho and phi are thus always the same, those for z range from 1 to the largest considered scale.
#' @examples
#' dta <- dt_analysis( blossom )
#' plot( dta )
#' dtb <- dt_analysis( blossom, return_fields=FALSE )
#' X11()
#' plot( dtb )
#' @seealso \code{\link{fld2dt}}, \code{\link{dt2cen}}
#' @export
dt_analysis <- function( fld, neg.rm=TRUE, rbr=NULL, pbr=NULL, zbr=NULL, mask=NULL, return_fields=TRUE, ... ){
    
    dt  <- fld2dt( fld, ... )
    if( neg.rm ) dt[ dt<0 ] <- 0
    
    cen <- dt2cen( dt, mask=mask )
    ms  <- dtmean( dt )
    ce  <- c( dt2cen( ms ) )
    names( ce ) <- c( "rho", "phi", "z" )
    res <- list( ms=ms, ce=ce )
    
    if( is.null(rbr) ) rbr <- 33
    if( is.null(pbr) ) pbr <- 33
    if( is.null(zbr) ) zbr <- 33
    
    if( length(rbr) == 1 ) rbr <- seq( 0,1,,rbr )
    if( length(pbr) == 1 ) pbr <- seq( 0,180,,pbr )
    if( length(zbr) == 1 ) zbr <- seq( 0,nrow(dt),,zbr )
    
    res$rhist <- hist( cen[ ,,1 ], breaks=rbr, plot=FALSE )$density
    res$phist <- hist( cen[ ,,2 ], breaks=pbr, plot=FALSE )$density
    res$zhist <- hist( cen[ ,,3 ], breaks=zbr, plot=FALSE )$density
    
    res$rmi <- ( rbr[-1] + rbr[ -length(rbr) ] )/2
    res$pmi <- ( pbr[-1] + pbr[ -length(pbr) ] )/2
    res$zmi <- ( zbr[-1] + zbr[ -length(zbr) ] )/2
    
    if( return_fields ){ 
        res$cen <- cen
        res$dt  <- dt
        res$uv  <- cen2uv(cen)
        res$fld <- fld
    }
    
    class( res ) <- "dtana"
    return( res )
}

