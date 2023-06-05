#' Colorimetric assay calibration
#'
#' Generates a standard curve for colormetric assays to determine macromolecular concentration.
#'
#'@param conc The concentration of the standards.
#'@param abs The measured absorbance values of the standards. Must be in the same order as the concentrations. Can be a concatenation of replicates if each replicate sequence is in the same order as the enumerated concentrations.
#'@param with.blank Boolean. Whether or not a blank (buffer only) is included in the values.
#'@param nrep The number of replicates included.
#'
#' @examples calibrate_colorimetric( conc, abs, with.blank == TRUE, nrep = 1 )
#'
#' @export
calibrate_colorimetric <- function( conc, abs, with.blank = TRUE, nrep = 1 ) {

  conc <- rep( conc, nrep )

  if ( isTRUE( with.blank ) ){
    abs <- matrix( abs,  ncol = nrep )
    blank <- abs[ which( rowSums( matrix( conc, ncol = nrep ) ) == 0 ), ]
    abs <- sweep( x = abs, MARGIN = 2, STATS = blank, FUN = "-" )
    abs <- abs[ which( rowSums( abs ) != 0 ), ]
    abs <- as.vector( abs )
    blank <- mean( blank )
  }

  if ( isFALSE( with.blank ) ){
    blank <- 0
  }

  conc <- conc[ conc != 0 ]
  linfit <- lm( abs ~ conc )

  coefs <- coefficients( linfit )
  coefs <- c( coefs, blank )
  names( coefs ) <- c( "intercept", "slope", "blank" )

  fit <- function( a, df = 1, return.conc = TRUE ) {
    if ( isTRUE( return.conc ) ){
      ( ( ( a - coefs[3] ) - coefs[1] ) / coefs[2] ) * df
      } else {
      coefs[1] + coefs[2] * a
        }
    }
  attr( fit, "coefficients" ) <- coefs

  return( fit )
}
