#' The 2D forward and inverse dualtree complex wavelet transform
#'
#' These functions perform the dualtree complex wavelet analysis and synthesis, either with or without decimation.
#' @param fld real matrix representing the field to be transformed
#' @param pyr a list containing arrays of complex coefficients for each level of the decomposition, produced by \code{dtcwt( ..., dec=TRUE )}
#' @param fb1 A list of filter coefficients for the first level. Currently only \code{near_sym_b} and \code{near_sym_b_bp} are implemented
#' @param fb2 A list of filter coefficients for all following levels. Currently only \code{qshift_b} and \code{qshift_b_bp} are implemented
#' @param J  number of levels for the decomposition. Defaults to \code{log2( min(Nx,Ny) )} in the decimated case and \code{log2( min(Nx,Ny) ) - 3} otherwise
#' @param dec whether or not the decimated transform is desired
#' @param verbose if TRUE, the function tells you which level it is working on
#' @return if dec=TRUE a list of complex coefficient fields, otherwise a complex \code{J x Nx x Ny x 6} array.
#' 
#' @details This is the 2D complex dualtree wavelet transform as described by Selesnick et al. (2005). It consists of four discerete wavelet transform trees, generated from two filter banks a and b by applying one set of filters to the rows and another (or the same) one to the columns. The 12 resulting coefficients are combined into six complex values representing six directions (15°, 45°, 75°, 105°, 135°, 165°). 
#' In the decimated case (dec=TRUE), each convolution is followed by a downsampling by two, meaining that the size of the six coefficient fields is cut in half at each level. The decimated transform can be reversed to recover the original image. For the \code{near_sym_b} and \code{qshift_b} filter banks, this reconstrcution should be basically perfect. In the case of the the \code{b_bp} filters, non-negligible artifacts appear near +-45° edges. 
#'
#' @note At present, the inverse transform only works if the input image had dimensions \code{2^N x 2^N}. You can use \code{\link{boundaries}} to achieve that.
#' 
#' @seealso \code{\link{filterbanks}}, \code{\link{fld2dt}}
#'
#' @references Kingsbury, Nick (1999) <doi:10.1098/rsta.1999.0447>.
#' Selesnick, I.W., R.G. Baraniuk, and N.C. Kingsbury (2005) <doi:10.1109/MSP.2005.1550194>
#' @author Nick Kingsbury (canonical MATLAB implementation), Rich Wareham (open source Python implementation, \url{https://github.com/rjw57/dtcwt}), Sebastian Buschow (R port).
#' @examples
#' oldpar <- par( no.readonly=TRUE )
#' # forward transform
#' dt <- dtcwt( blossom )
#' par( mfrow=c(2,3), mar=rep(2,4) )
#' for( j in 1:6 ){
#'     image( blossom, col=grey.colors(32,0,1) )
#'     contour( Mod( dt[[3]][ ,,j ] )**2, add=TRUE, col="green" )
#' } 
#' par( oldpar ) 
#'
#' # exmaple for the inverse transform
#' blossom_i <- idtcwt( dt )
#' image( blossom - blossom_i )
#' 
#' # example for a non-square case
#' boy <- blossom[50:120, 50:150]
#' bc  <- put_in_mirror(boy, 128)
#' dt  <- dtcwt(bc$res)
#' idt <- idtcwt(dt)[ bc$px, bc$py ]
#' @name dualtree-transform
NULL

#' @export
#' @rdname dualtree-transform
dtcwt <- function( fld, fb1=near_sym_b, fb2=qshift_b, J=NULL, dec=TRUE, verbose=FALSE ){


    mc <- my_conv

    # check whether either of the filter banks has additional filters for the H-Hi component
    sym1 <- length( fb1 ) > 4 
    sym2 <- length( fb2 ) > 8 
    
    # prepare dimensions, make sure we have an appropriate rectangle
    Nx <- nrow( fld )
    Ny <- ncol( fld )
    
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
        pass01 <- mc( fld, get(fil1), dec=FALSE )
        if( sp1 ) pass01b <- mc( fld, Hi2, dec=FALSE )
        
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



#' @export
#' @rdname dualtree-transform
idtcwt <- function( pyr, fb1=near_sym_b, fb2=qshift_b, verbose=TRUE ){

    nx <- nrow( pyr[[1]] )
    ny <- ncol( pyr[[1]] )

    if( abs( floor( log2( nx ) ) - log2(nx) ) > 1e-10 | nx!=ny ) stop( "The inverse transform currently only works for squares of size 2^N x 2^N. Sorry." )

    mc <- function( x,y ) return( my_conv( upsample( x ), y, dec=FALSE, inverse=TRUE ) )

    mc <- list( a=mc, b=mc )
    
    sym1 <- length( fb1 ) > 4 
    sym2 <- length( fb2 ) > 8
    J    <- length( pyr ) - 1
    
    # get upper level filters
    Lo  <- list( a=fb2$g0a, b=fb2$g0b )
    Hi  <- list( a=fb2$g1a, b=fb2$g1b )
    Hi2 <- list( a=fb2$g2a, b=fb2$g2b )

    res <- list( LoaLob = as.matrix( pyr[[J+1]][1] ),
                 LobLoa = as.matrix( pyr[[J+1]][2] ),
                 LoaLoa = as.matrix( pyr[[J+1]][3] ),
                 LobLob = as.matrix( pyr[[J+1]][4] ) 
               )
    
    # loop over all levels
    for( j in J:1 ){
        tmp <- c2q( pyr[[j]] )
        
        if( j==1 ){
            # get level 1 filters when we need them ... 
            Lo   <- list( a=c( fb1$g0o, 0 ), b=c( fb1$g0o, 0 ) )
            Hi   <- list( a=fb1$g1o, b=fb1$g1o )
            Hi2  <- list( a=fb1$g2o, b=fb1$g2o )
            mc$b <- function( x, y ) return( my_conv( shift1( upsample( x ), forward=FALSE ),y, dec=FALSE, inverse=TRUE ) )
        }
        
        # loop over all combinations (aa,ab,ba,bb) to acces all four trees
        for( f1 in c("a","b") ){
            for( f2 in c("a","b") ){
            
                # get the names of the four branches
                LoLo <- paste0( "Lo",f1,"Lo",f2 )
                HiLo <- paste0( "Hi",f1,"Lo",f2 )
                LoHi <- paste0( "Lo",f1,"Hi",f2 )
                HiHi <- paste0( "Hi",f1,"Hi",f2 )
                
                # undo the second level
                y1 <- t( mc[[f2]]( t( res[[LoLo]] ),  Lo[[f2]]  ) + 
                         mc[[f2]]( t( tmp[[LoHi]] ),  Hi[[f2]]  ) ) 
                
                # undo the first level
                if( (j==1 & sym1)  | (j>1 & sym2) ){
                    y2  <- mc[[f2]]( t( tmp[[HiLo]] ),  Lo[[f2]]  )
                    y2b <- mc[[f2]]( t( tmp[[HiHi]] ),  Hi2[[f2]]  )
                    
                    y   <- mc[[f1]]( y1, Lo[[f1]] ) + 
                           mc[[f1]]( t(y2), Hi[[f1]] ) + 
                           mc[[f1]]( t(y2b), Hi2[[f1]] )
                    
                }else{         
                    y2 <- t( mc[[f2]]( t( tmp[[HiLo]] ),  Lo[[f2]]  ) + 
                             mc[[f2]]( t( tmp[[HiHi]] ),  Hi[[f2]]  ) ) 
                            
                    y <- mc[[f1]]( y1, Lo[[f1]] ) + mc[[f1]]( y2, Hi[[f1]] )
                }
                res[[LoLo]] <- y
            }
        }
    }
    res <- res[[1]] + res[[2]] + res[[3]] + res[[4]]
    return( Re(res) )
}
