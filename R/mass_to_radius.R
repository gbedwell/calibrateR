#' Convert molecular weight values to hydrodynamic radius
#'
#'
#'@param mass The void volume of the column.
#'
#'
#' @examples calibrate_sec( masses, void, cv, vols )
#'
#' @export
mass_to_radius <- function( masses ) {
  rh <- ( 10^( ( -0.204 + ( 0.357 * log10( masses ) ) ) ) ) / 10
  return( rh )
}
