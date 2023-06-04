#' Calculate Gibson assembly reaction volumes
#'
#'@param insert.conc Vector of insert concentrations.
#'@param insert.len Vector of insert lengths in bp.
#'@param vec.conc Vector concentration. Must be a single value.
#'@param vec.len Vector length. Must be a single value.
#'@param molar.ratio Molar ratio of insert:vector
#'@param vec.mass Mass of vector
#'@param final.vol Final volume of the reaction before addition of 2x master mix.
#'
#'
#'@examples gibson()
#'
#'@export
gibson <- function( insert.conc, insert.len, vec.conc, vec.len,
                    molar.ratio = 3, vec.mass = 50, final.vol = 10 ) {

  if ( length( vec.len ) > 1 | length( vec.conc ) > 1 ){
    stop( "Only a single vector input is allowed.", call. = FALSE )
  }

  if ( length( insert.conc ) != length( insert.len ) ){
    stop( "Insert concentration and insert length vectors contain different numbers of terms.", call. = FALSE )
  }

  nfrag <- length( insert.conc )

  pmol.vec <- ( vec.mass * 1000 ) / ( vec.len * 650 )
  vol.vec <- round( vec.mass / vec.conc, 2 )

  pmol.insert <- pmol.vec * molar.ratio
  mass.insert <- ( pmol.insert *  ( insert.len * 650 ) ) / 1000

  vol.insert <- round( mass.insert / insert.conc, 2 )
  tot.insert.vol <- sum( vol.insert )
  names( vol.insert ) <- paste( "Fragment", 1:nfrag, sep = " " )

  vol.water <- round( final.vol - ( vol.vec + tot.insert.vol ), 2 )

  vol.mix <- final.vol

  df <- data.frame( Volume = vol.insert )

  df <- rbind( df,
               data.frame(
                 Volume = c( Vector = vol.vec,
                             Water = vol.water,
                             `Master Mix` = vol.mix ) ) )

  return( df )

}
