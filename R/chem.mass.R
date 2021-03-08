#' Calculates the mass any chemical to add to a solution to reach a target concentration
#'
#'
#'
#'@param chemicals The name of each chemical
#'@param mass The mass of the chemical
#'@param final.conc The final concentration of the solution in MOLAR units
#'@param final.vol The final volume of solution
#'@param vol.unit The desired units of output
#'
#'
#' @examples chem.mass()
#'
#' @export
chem.mass <- function(chemicals, mass, final.conc, final.vol, vol.unit = c("mL","L"), conc.unit = c("M","mM")){

  dat <- data.frame(chemical = chemicals)

  if(vol.unit == "mL"){

    if(conc.unit == "M"){
      dat <- dat %>%
        dplyr::mutate(g.to.add = mass*final.conc*(final.vol/1000),
                      final.volume = final.vol) }

    if(conc.unit == "mM"){
      dat <- dat %>%
        dplyr::mutate(g.to.add = mass*(final.conc/1000)*(final.vol/1000),
                      final.volume = final.vol) }


  if(vol.unit == "L"){
    dat <- dat %>%
      dplyr::mutate(g.to.add = mass*final.conc*final.vol,
                    final.vol = final.vol*1000) }

    if(conc.unit == "M"){
      dat <- dat %>%
        dplyr::mutate(g.to.add = mass*final.conc*final.vol,
                      final.volume = final.vol) }

    if(conc.unit == "mM"){
      dat <- dat %>%
        dplyr::mutate(g.to.add = mass*(final.conc/1000)*final.vol,
                      final.volume = final.vol) }


  return(dat) } }
