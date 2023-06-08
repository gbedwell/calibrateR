#' Calculate CFU/mL dilutions
#'
#'@param od The measured OD of the starter culture.
#'@param od.df The dilution factor of the OD measurement.
#'@param count.cf The correction factor that relates OD to CFU/mL. Defaults to 8E8.
#'@param target.cpmL The desired CFU/mL in the new culture.
#'@param target.vol The desired volume of the new culture.
#'
#'
#' @examples cfu()
#'
#' @export
cfu <- function( od, od.df, count.cf = 8E8, target.cpmL, target.vol ){

  corr.od <- od * od.df
  cpmL <- corr.od * count.cf
  vol <- ( target.vol * target.cpmL ) / cpmL
  return( vol )

}
