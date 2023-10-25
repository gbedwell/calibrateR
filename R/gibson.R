#' Calculate Gibson assembly reaction volumes
#'
#'@param insert.conc Vector of insert concentrations.
#'@param insert.len Vector of insert lengths in bp.
#'@param vec.conc Vector concentration. Must be a single value.
#'@param vec.len Vector length. Must be a single value.
#'@param n.frags The number of fragments in each reaction.
#'@param molar.ratio Molar ratio of insert:vector
#'@param vec.mass Mass of vector
#'@param final.vol Final volume of the reaction before addition of 2x master mix.
#'@param ids Reaction identifiers. Can be NULL.
#'
#'
#'@examples gibson()
#'
#'@export
gibson <- function( insert.conc, insert.len, vec.conc, vec.len, n.frags,
                    molar.ratio = 3, vec.mass = 50, final.vol = 5, ids = NULL ) {

  if( length( vec.conc ) == 1 && length( vec.len ) == 1 &&
      length( molar.ratio ) == 1 && length( vec.mass ) == 1 &&
      length( final.vol ) == 1 && length( n.frags ) != 1 ){

    vec.conc <- rep( vec.conc, length( n.frags ) )
    vec.len <- rep( vec.len, length( n.frags ) )
    molar.ratio <- rep( molar.ratio, length( n.frags ) )
    vec.mass <- rep( vec.mass, length( n.frags ) )
    final.vol <- rep( final.vol, length( n.frags ) )

    }

  if ( length( unique( c( length( vec.conc ), length( vec.len ),
                          length( molar.ratio ),length( vec.mass ),
                          length( final.vol ), length( n.frags ) ) ) ) != 1 ){

    stop( "vec.conc, vec.len, n.frags, molar.ratio, vec.mass, and final.vol vectors must be the same length.",
          call. = FALSE)

      }

  if ( !is.null( ids ) ) {
    if ( length( ids ) != length( n.frags ) ){
      stop( "ids vector must have the same length as n.frags vector.", call. = FALSE )
    }
  }

  if ( length( insert.conc ) != length( insert.len ) ){
    stop( "insert.conc and insert.len vectors must be the same length.", call. = FALSE )
  }

  x <- 0

  ll <- list()

  for ( i in seq_along( n.frags ) ){
    f <- n.frags[i]
    vc <- vec.conc[i]
    vl <- vec.len[i]
    mr <- molar.ratio[i]
    vm <- vec.mass[i]
    fv <- final.vol[i]

    start <- x + 1
    end <- x + f

    ic <- insert.conc[ start:end ]
    il <- insert.len[ start:end ]

    pmol.vec <- ( vm * 1000 ) / ( vl * 650 )
    vol.vec <- round( vm / vc, 2 )

    pmol.insert <- pmol.vec * mr
    mass.insert <- ( pmol.insert *  ( il * 650 ) ) / 1000

    vol.insert <- round( mass.insert / ic, 2 )
    tot.insert.vol <- sum( vol.insert )
    names( vol.insert ) <- paste( "Fragment", 1:f, sep = " " )

    vol.water <- round( fv - ( vol.vec + tot.insert.vol ), 2 )

    vol.mix <- fv

    df <- data.frame( Volume = vol.insert )

    df <- rbind( df,
                 data.frame(
                   Volume = c( Vector = vol.vec,
                               Water = vol.water,
                               `Master Mix` = vol.mix ) ) )

    x <- x + f

    ll <- c( ll, list( df ) )
  }

  if ( !is.null( ids ) ){
    names( ll ) <- ids
  }

  return( ll )

}
