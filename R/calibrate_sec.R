#' Gel filtration column calibration
#'
#' Calculates the calibration curve for a size exclusion column using standards of known mass. Standard mass values should be entered in kDa. Hydrodynamic radius values are calculated from the mass values in nm units. The hydrodynamic radius - molar mass relationship is defined in https://doi.org/10.1046/j.0014-2956.2001.02649.x
#'
#'@param vols The elution volumes of the standards.
#'@param masses The known masses of each of the standards.
#'@param parameter One of "mw" or "rh". Defines whether the column is being calibrated for molecular weight or hydrodynamic radius.
#'@param normalized Boolean. Whether or not the input elution volumes should be normalized relative to the given void and column volumes.
#'@param void The void volume of the column.
#'@param cv The column volume.
#'
#'@importFrom dplyr arrange
#'@importFrom dplyr mutate
#'@import nls.multstart
#'
#' @examples calibrate_sec( masses, void, cv, vols )
#'
#' @export
calibrate_sec <- function( vols, masses, parameter = c( "mw", "rh" ),
                           void = NULL, cv = NULL, normalized = FALSE ) {

  if (isTRUE( normalized ) ){
    if ( is.null( void ) | is.null( cv ) ){
      stop( "void and cv must be defined to normalize elution volumes.", call. = FALSE )
    }
  }

  parameter <- match.arg( parameter )

  if ( parameter != "mw" & parameter != "rh" ){
    stop( "parameter must be defined as either 'mw' or 'rh'." )
  }

  df <- data.frame( vol = vols,
                    mw = masses ) |>
    dplyr::mutate(rh = mass_to_radius( masses = mw ) )

  if ( isTRUE( normalized ) ){
    df <- df |>
      dplyr::mutate( vol = normalize_ev( vols = vol, void = void, cv = cv ) )
  }

  if ( parameter == "mw" ){
    fit <- nls.multstart::nls_multstart(log( mw ) ~ a * vol^3 + b * vol^2 + c * vol^1 + d,
                                        data = df,
                                        iter = 250,
                                        start_lower = c(a = -250, b = -250, c = -250, d = -250),
                                        start_upper = c(a = 250, b = 250, c = 250, d = 250) )
  }

  if ( parameter == "rh" ){
    fit <- nls.multstart::nls_multstart(log( rh ) ~ a * vol^3 + b * vol^2 + c * vol^1 + d,
                                        data = df,
                                        iter = 250,
                                        start_lower = c(a = -250, b = -250, c = -250, d = -250),
                                        start_upper = c(a = 250, b = 250, c = 250, d = 250) )
  }

  coefs <- coefficients( fit )

  trend <- function( v ) { exp( coefs[1] * v^3 + coefs[2] * v^2 + coefs[3] * v + coefs[4] ) }
  attr( trend, "coefficients" ) <- coefs

  return( trend )
}

