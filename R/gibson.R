#' Calculates the volumes of vector and insert of a given series of Gibson assembly reactions
#'
#'
#'
#'@param insert.list List of all inserts being put into plasmids
#'@param vector.list List of all vectors being cloned into. Must enter one vector for each insert, even if there are duplicates.
#'@param insert.conc The concentration of insert in ng/µL
#'@param insert.len The insert size in bp
#'@param vec.conc The concentration of the vector in ng/µL
#'@param vec.len The linearized vector size in bp
#'@param molar.ratio The molar ratio of insert to vector
#'@param final.vol The final volume of the reaction BEFORE addition of 2x reaction buffer.
#'@param vec.mass The target mass of vector. If none given, assumed to be 50 ng.
#'@param combine.fragments \code{TRUE} or \code{FALSE}. Tells the function whether or not to consider all fragments as going into a single reaction.
#'
#'
#' @examples gibson()
#'
#' @export
gibson <- function(insert.list, vector.list, insert.conc, insert.len, vec.conc, vec.len, molar.ratio, final.vol, vec.mass=50, combine.fragments=FALSE) {

  dat <- data.frame(insert = insert.list,
                    vector = vector.list)

  if(!isTRUE(combine.fragments)){
    dat <- dat %>%
      dplyr::mutate(pmol.vec = (vec.mass*1000)/(vec.len*650),
                    vol.vec = round(50/vec.conc, 2),
                    pmol.insert = pmol.vec*molar.ratio,
                    mass.insert = (pmol.insert*(insert.len*650))/1000,
                    vol.insert = round(mass.insert/insert.conc, 2),
                    vol.water = round(final.vol - (vol.vec+vol.insert), 2),
                    vol.mix = final.vol) %>%
      dplyr::select(insert, vector, vol.vec, vol.insert, vol.water, vol.mix) %>%
      magrittr::set_colnames(c("Insert","Vector","Vector Volume (µL)","Insert Volume (µL)",
                               "Volume Water (µL)","Volume MM (µL)"))
    }
  else{
    dat <- dat %>%
      dplyr::mutate(pmol.vec = (vec.mass*1000)/(vec.len*650),
                    vol.vec = round(50/vec.conc, 2),
                    pmol.insert = pmol.vec*molar.ratio,
                    mass.insert = (pmol.insert*(insert.len*650))/1000,
                    vol.insert = round(mass.insert/insert.conc, 2),
                    vol.water = round(final.vol - (vol.vec+sum(vol.insert)), 2),
                    vol.mix = final.vol) %>%
      dplyr::select(insert, vector, vol.vec, vol.insert, vol.water, vol.mix) %>%
      magrittr::set_colnames(c("Insert","Vector","Vector Volume (µL)","Insert Volume (µL)",
                               "Volume Water (µL)","Volume MM (µL)"))
    }
  return(dat)
}
