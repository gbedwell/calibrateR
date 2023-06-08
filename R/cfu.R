#' Calculate CFU/mL dilutions
#'
#'@param od The measured OD of the starter culture.
#'@param od.df The dilution factor of the OD measurement.
#'@param count.cf The correction factor that relates OD to CFU/mL. Defaults to 8E8.
#'@param target.cfu The desired CFU/mL in the new culture.
#'@param target.vol The desired volume of the new culture.
#'
#'
#' @examples cfu()
#'
#' @export
cfu <- function( od, od.df, count.cf = 8E8, target.cfu, target.vol ){

  if ( length( od.df > 1 ) & length( od.df ) != length( od ) ){
    stop("It looks like you have multiple dilution factors, but the number of dilution factors does not match the number of OD readings!", call. = FALSE)
  }

  corr.od <- od * od.df
  cfu <- corr.od * count.cf
  vol <- ( target.vol * target.cfu ) / cfu
  return( vol )

}
