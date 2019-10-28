#' The 2D inverse dualtree complex wavelet transform
#'
#' Reconstructs an image from the pyramid of complex directional wavelet coefficients.
#'
#' @param pyr a list containing arrays of complex coefficients for each level of the decomposition, produced by \code{dtcwt( ..., dec=TRUE )}.
#' @param fb1 the filter bank for the first level
#' @param fb2 the filter bank for all following levels
#' @param verbose if true, the function will say a few words while doing its thing. 
#' @param boundaries how to handle the boundary conditions, should be the same as for the decomposition.
#' @return a real array of size \code{ 2N x 2M } where \code{ dim( pyr[[1]] ) = (M,N,6) }.
#' @details This function re-arranges the six complex daughter coefficients back into the four trees, convolves them with the synthesis wavelets and adds everything up to recover an image. For the \code{near_sym_b} and \code{qshift_b} filter banks, this reconstrcution should be basically perfect. In the case of the the \code{b_bp} filters, non-negligible artifacts appear near +-45° edges. 
#' @note At present, only boundaries="periodic" actually works :(
#' @seealso \code{\link{dtcwt}}
#'
#' @references Selesnick, I.W., R.G. Baraniuk, and N.C. Kingsbury. “The Dual-Tree Complex Wavelet Transform.” IEEE Signal Processing Magazine 22, no. 6 (November 2005): 123–51. \url{https://doi.org/10.1109/MSP.2005.1550194}.
#' @examples
#' py <- dtcwt( boys )
#' boys_i <- idtcwt( py )
#' image( boys - boys_i )
#' @export
idtcwt <- function( pyr, fb1=near_sym_b, fb2=qshift_b, verbose=TRUE, boundaries="periodic" ){

    mc <- function( x,y ) return( my_conv( upsample( x ), y, dec=FALSE, boundaries=boundaries ) )

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
            mc$b <- function( x, y ) return( my_conv( shift1( upsample( x ), forward=FALSE ),y, dec=FALSE, boundaries=boundaries ) )
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
    return( res )
}
