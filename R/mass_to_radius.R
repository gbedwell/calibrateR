#' Convert molecular weight values to hydrodynamic radius
#'
#'
#'@param mass The void volume of the column.
#'
#'
#' @examples mass_to_radius( masses = c( 1E4, 3.7E4 ) )
#'
#' @export
mass_to_radius <- function( masses ) {
  rh <- ( 10^( ( -0.204 + ( 0.357 * log10( masses ) ) ) ) ) / 10
  return( rh )
}
