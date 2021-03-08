#' Calculates the volume of each component in a dilution
#'
#'
#'
#'@param samples The name of each sample being diluted
#'@param final.vol The targeted final volume
#'@param final.conc The targeted final concentration
#'@param starting.conc The starting concentration
#'@param combine Whether or not the final mixture should be a mixture of two analytes (e.g. for a binding reaction)
#'
#'
#' @examples dilution.calc()
#'
#' @export
dilution.calc <- function(samples, final.vol, final.conc, starting.conc, combine=c("yes","no")){

  if(combine == "no"){

    dat <- data.frame(species = samples)

    dat <- dat %>%
      dplyr::mutate(vol.sample = (final.vol*final.conc)/starting.conc,
                    vol.buffer = final.vol - vol.sample) }

  if(combine == "yes"){

    dat <- data.frame(species = samples)

    dat <- dat %>%
      dplyr::mutate(vol.sample = (final.vol*final.conc)/starting.conc,
                    vol.buffer = final.vol - sum(vol.sample)) }

  return(dat) }
