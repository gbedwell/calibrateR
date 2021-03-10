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
#'
#'
#' @examples gibson()
#'
#' @export
gibson <- function(insert.list, vector.list, insert.conc, insert.len, vec.conc, vec.len, molar.ratio, final.vol, vec.mass=NULL){

  if(is.null(vec.mass)){

    writeLines(c("--------------------------------------------------------",
                 "Assuming the desired vector mass is 50 ng",
                 "--------------------------------------------------------",
                 ""))

    dat <- data.frame(sample = insert.list,
                      vector = vector.list)

    dat <- dat %>%
      dplyr::mutate(pmol.vec = (50*1000)/(vec.len*650),
                    vol.vec = 50/vec.conc,
                    pmol.insert = pmol.vec*molar.ratio,
                    mass.insert = (pmol.insert*(insert.len*650))/1000,
                    vol.insert = mass.insert/insert.conc,
                    vol.water = final.vol - (vol.vec+vol.insert),
                    vol.mix = final.vol) %>%
      dplyr::select(sample, vector, vol.vec, vol.insert, vol.water, vol.mix)

    return(dat) }

  else{
    dat <- data.frame(sample = insert.list,
                      vector = vector.list)

    dat <- dat %>%
      dplyr::mutate(pmol.vec = (vec.mass*1000)/(vec.len*650),
                    vol.vec = vec.mass/vec.conc,
                    pmol.insert = pmol.vec*molar.ratio,
                    mass.insert = (pmol.insert*(insert.len*650))/1000,
                    vol.insert = mass.insert/insert.conc,
                    vol.water = final.vol - (vol.vec+vol.insert),
                    vol.mix = final.vol) %>%
      dplyr::select(sample, vector, vol.vec, vol.insert, vol.water, vol.mix)

    return(dat) } }
