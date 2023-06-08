#' Normalize SEC elution volumes
#'
#' Normalize SEC elution volumes to be between 0 and 1 using measured void volume and column volume
#'
#'@param void The void volume of the column.
#'@param cv The column volume.
#'
#' @examples calibrate_sec( masses, void, cv, vols )
#'
#' @export
normalize_ev <- function( vols, void, cv ) {
  vols = ( vols - void ) / ( cv - void )
  return( vols )
}
