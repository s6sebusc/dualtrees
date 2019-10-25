#' Transform data from the four trees to six fields of complex coefficients.
#'
#' This function takes the four combinations of high- and low passes from the four trees (ab, ba, aa, bb) and re-arranges them into the six directional complex daughter wavelets.
#' @param q a list of wavelet coefficients named LoaHia, LobHib, HiaLoa, ...
#' @return a complex array of size nx, ny, 6
#' @keywords drudenfuss
#' @export
#' @examples
#' q2c( q )

q2c <- function( q ){
    comp <- with( q,{
        res <- array( dim=c( dim( LoaHia ), 6 ) , data=NA )
        res[ ,,1 ] <- LoaHia - LobHib + 1i*( LobHia + LoaHib )
        res[ ,,6 ] <- LoaHia + LobHib + 1i*( LobHia - LoaHib )
        res[ ,,3 ] <- HiaLoa - HibLob + 1i*( HibLoa + HiaLob )
        res[ ,,4 ] <- HiaLoa + HibLob + 1i*( HibLoa - HiaLob )
        res[ ,,5 ] <- HiaHia - HibHib + 1i*( HibHia + HiaHib )
        res[ ,,2 ] <- HiaHia + HibHib + 1i*( HibHia - HiaHib )
        res / sqrt(2) 
    } )
    return( comp )
}

#' Transform six fields of complex coefficients back into four trees.
#'
#' This function takes the the six directional complex daughter wavelet coefficients and re-constructs the three combinations of high- and low passes from the four trees (ab, ba, aa, bb).
#' @param comp complex array of dimnesions nx, ny, 6
#' @return a list of low- and high-pass components from the four trees, names LoaHia, LobHib, etc.
#' @keywords drudenfuss
#' @export
#' @examples
#' c2q( comp )
c2q <- function( comp ){
    ra <- Re( comp ) / sqrt(2)
    rb <- Im( comp ) / sqrt(2)
    q  <- list( LoaHia = ra[,,6] + ra[,,1],
                LobHib = ra[,,6] - ra[,,1], 
                HiaLoa = ra[,,4] + ra[,,3],
                HibLob = ra[,,4] - ra[,,3], 
                HiaHia = ra[,,2] + ra[,,5],
                HibHib = ra[,,2] - ra[,,5],
                
                LobHia = rb[,,1] + rb[,,6],
                LoaHib = rb[,,1] - rb[,,6], 
                HibLoa = rb[,,3] + rb[,,4],
                HiaLob = rb[,,3] - rb[,,4], 
                HibHia = rb[,,5] + rb[,,2],
                HiaHib = rb[,,5] - rb[,,2]
    )
    return( q )
}

#' delete every second row of a matrix
#' @export
decimate <- function( mat, odd=FALSE, dec=TRUE ){
    # if dec = FALSE, this does straight nothing :D
    # odd =TRUE/FALSE determines which samples are kept
    if( dec ){
        if( !odd ){ #!
            mat <- mat[ seq( 1, nrow(mat), 2 ), ,drop=FALSE ]
        }else{
            mat <- mat[ seq( 2, nrow(mat), 2 ), ,drop=FALSE ]
        }
    }
    return( mat )
}

#' add rows with zeroes to a matrix
#' @export
upsample <- function( mat, odd=TRUE ){
    mat <- as.matrix( mat )
    res <- matrix( nrow=2*nrow(mat), ncol=ncol(mat), data=0 )
    res[ seq( if(odd) 2 else 1, nrow(res),2 ), ] <- mat  
    return( res )
}

#' shift a matrix forward or backward by one row
#' @export
shift1 <- function( x, forward=TRUE ){
    x <- as.matrix( x )
    if( forward ){
        x <- rbind( x[ nrow(x),,drop=FALSE ], x[ -nrow(x),,drop=FALSE ] )
    }else{
        x <- rbind( x[-1,,drop=FALSE], x[1,,drop=FALSE] )
    }
#     
    return( x )
}


#' insert holes into a filter?
#' @export
holes <- function( fil, second=TRUE ){
    res <- rbind( fil, rep( 0, length(fil) ) ) 
    # add a zero at the other end to keep the wavelet symmetrical!
    return( c( 0, res ) )
}
