#' Determine macromolecule concentration from UV/Vis measurements.
#'
#'@param abs The measured absorbance value at the wavelength of interest.
#'@param ext The extinction coefficient of the molecule. The units of the extinction coefficient determine the units of concentration.
#'@param df The dilution factor of the measurement.
#'@param plen The path length of the measurement. Defaults to 1 cm.
#'@param type The type of molecule being analyzed.
#'
#'
#' @examples uv_vis()
#'
#' @export
uv_vis <- function( abs, ext = NULL, df = 1, plen = 1,
                    type = c("protein","dsDNA","ssDNA","ssRNA") ){

  type <- match.arg( type )

  if ( type != "protein" & type != "dsDNA" & type != "ssDNA" & type != "ssRNA" ){
    stop( "Type must be one of 'protein', 'dsDNA', 'ssDNA', or 'ssRNA'", call. = FALSE )
  }

  if ( is.null( ext ) ){
    if ( type == "protein" ){
      stop( "Extinction coefficient must be defined for proteins.", call. = FALSE )
      } else if ( type == "dsDNA" ){
        ext <- 50
        } else if ( type == "ssDNA" ){
        ext <- 33
        } else if ( type == "ssRNA" ){
        ext <- 40
        }
  }

  conc <- ( abs / ( ext * plen ) ) * df

  return( conc )

}
