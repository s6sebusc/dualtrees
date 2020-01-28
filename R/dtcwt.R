#' The 2D forward dualtree complex wavelet transform
#'
#' This function performs the dualtree complex wavelet analysis, either with or withour decimation
#' @param mat the real matrix we wish to transform
#' @param fb1 A list of analysis filter coefficients for the first level. Currently only near_sym_b and near_sym_b_bp are implemented
#' @param fb2 A list of analysis filter coefficients for all following levels. Currently only qshift_b and qshift_b_bp are implemented
#' @param J  number of levels for the decomposition. Defaults to \code{log2( min(Nx,Ny) )} in the decimated case and \code{log2( min(Nx,Ny) ) - 3} otherwise
#' @param dec whether or not the decimated transform is desired
#' @param verbose if TRUE, the function tells you which level it is working on
#' @return if dec=TRUE a list of complex coefficient fields, otherwise a complex J * Nx * Ny * 6 array.
#' 
#' @details This is the 2D complex dualtree wavelet transform as described by Selesnick et al 2005. It consists of four discerete wavelet transform trees, generated from two filter banks a and b by applying one set of filters to the rows and the same ot the other to the columns. 
#' In the decimated case (dec=TRUE), each convolution is followed by a downsampling, meaining that the size of the six coefficient fields is cut in half at each level. 
#'
#' @note The inverse transform only works if the input is of size \code{2^N x 2^N}.
#' 
#' @seealso \code{\link{idtcwt}}
#'
#' @references Selesnick, I.W., R.G. Baraniuk, and N.C. Kingsbury. “The Dual-Tree Complex Wavelet Transform.” IEEE Signal Processing Magazine 22, no. 6 (November 2005): 123–51. \url{https://doi.org/10.1109/MSP.2005.1550194}.
#' @examples
#' dt <- dtcwt( blossom )
#' par( mfrow=c(2,3), mar=rep(2,4) )
#' for( j in 1:6 ){
#'     image( blossom, col=grey.colors(32,0,1) )
#'     contour( Mod( dt[[3]][ ,,j ] )**2, add=TRUE, col="green" )
#' } 
#' @export
dtcwt <- function( mat, fb1=near_sym_b, fb2=qshift_b, J=NULL, dec=TRUE, verbose=TRUE ){


    mc <- my_conv

    # check whether either of the filter banks has additional filters for the H-Hi component
    sym1 <- length( fb1 ) > 4 
    sym2 <- length( fb2 ) > 8 
    
    # prepare dimensions, make sure we have an appropriate rectangle
    Nx <- nrow( mat )
    Ny <- ncol( mat )
    
    if( !dec ){
        if( max( abs( log2(Nx) - round( log2(Nx) ) ),
                 abs( log2(Ny) - round( log2(Ny) ) ) ) > 1e-6 ){ 
            stop( "input dimensions must be 2^N x 2^M" )
        }
    }
    
    # find the maximum level
    if( is.null( J ) ){
        if( dec ){
            J <- log2( min(Nx,Ny) )
        }else{
            J <- log2( min(Nx,Ny) ) - 3
        }
    }
    
    # get first level filters
    Hi <- c( fb1$h1o, 0 )  
    Lo <- fb1$h0o  
    if( sym1 ) Hi2 <- c( fb1$h2o , 0 )
    
    # perpare list for results
    res <- list()

    # start with level 1
    if( verbose ) print( 1 )
    fun <- list( a=identity, b=shift1 )
    
    tmp <- list()
    # loop over all combinations of filter a and b, high-pass and low-pass,
    # giving us four trees (aa,ab,ba,bb) with three daughters (hl,hh,lh) and 
    # their father (ll) who is passed to the next levels
    for( fil1 in c( "Lo", "Hi" ) ){
        sp1    <- fil1=="Hi" & sym1
        pass01 <- mc( mat, get(fil1), dec=FALSE )
        if( sp1 ) pass01b <- mc( mat, Hi2, dec=FALSE )
        
        for( f1 in c("a", "b") ){
            pass1 <- t( decimate( fun[[f1]]( pass01 ), dec=dec ) )
            if( sp1 ) pass1b <- t( decimate( fun[[f1]]( pass01b ), dec=dec ) )
            
            for( fil2 in c( "Lo", "Hi" ) ){
                sp2    <- fil2=="Hi" & sym1
                pass02 <- mc( pass1, get(fil2 ), dec=FALSE )
                if( sp1 & sp2 ) pass02b <- mc( pass1b, Hi2, dec=FALSE )
                
                for( f2 in c("a", "b") ){
                    if( sp1 & sp2 ){
                        pass2 <- t( decimate( fun[[f2]]( pass02b ), dec=dec ) )
                    }else{
                        pass2 <- t( decimate( fun[[f2]]( pass02 ), dec=dec ) )
                    }
                    tmp[[ paste0( fil1,f1,fil2,f2 ) ]] <- pass2
                }
            }
        }
    }
    res[[1]] <- q2c( tmp )
    
    # get upper level filters
    Lo  <- list( a=fb2$h0a, b=fb2$h0b )
    Hi  <- list( a=fb2$h1a, b=fb2$h1b )
    Hi2 <- list( a=fb2$h2a, b=fb2$h2b )
    
    for( j in 2:J ){
        if( verbose ) print( j )
        
        if( !dec ){ # in the undecimated case, add holes
            Lo  <- lapply( Lo, holes )
            Hi  <- lapply( Hi, holes )
            Hi2 <- lapply( Hi2, holes )
        }
        
        tmp2 <- list()
        # loop over all combinations of filter a and b, giving us again
        # four trees (aa,ab,ba,bb)
        for( f1 in c( "a", "b" ) ){
            for( f2 in c( "a", "b" ) ){
                LoLo  <- tmp[[ paste0( "Lo",f1,"Lo",f2 ) ]]
                
                # loop over all combinations of high- and low-pass, giving us
                # three daughters and one father, just as before
                for( fil1 in c( "Lo", "Hi" ) ){
                    sp1   <- fil1=="Hi" & sym2
                    pass1 <- t( mc( LoLo, get(fil1)[[f1]], dec=dec ) )
                    if( sp1 ) pass1b <- t( mc( LoLo, Hi2[[f1]], dec=dec ) )
                    
                    for( fil2 in c( "Lo", "Hi" ) ){
                        if( sp1 & fil2=="Hi" ){
                            pass2 <- t( mc( pass1b, Hi2[[f2]], dec=dec ) )
                        }else{
                            pass2 <- t( mc( pass1, get(fil2)[[f2]], dec=dec ) )
                        }
                        tmp2[[ paste0( fil1,f1,fil2,f2 ) ]] <-  pass2 
                    }
                }
            }
        }
        res[[j]] <- q2c( tmp2 )
        tmp <- tmp2
    }
    # save the final low pass
    res[[ J+1 ]] <- with( tmp, c( LoaLob, LobLoa, LoaLoa, LobLob ) )
    
    if( !dec ){ # re-save into an array if we did not decimate
        res_ar <- array( dim=c( J, Nx, Ny, 6 ), data=NA )
        for( j in 1:J ) res_ar[j,,,] <- res[[j]] 
        res <- res_ar
        class( res ) <- "dt"
    }else{
        class( res ) <- "dt_pyr"
    }
    
    return( res )
}
