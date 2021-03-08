#' Calculates the concentration of protein, dsDNA, ssDNA, or ssRNA from absorbance values at 280 or 260 nm
#'
#'
#'
#'@param type The type of molecule being analyzed \code{protein},\code{dsDNA},\code{ssDNA},or \code{ssRNA}
#'@param A260 The absorbance value at 260 nm
#'@param A280 The absorbance value at 280 nm
#'@param dil.fac The dilution factor of the measured sample
#'@param path.length The path length of the instrument in cm. Usually 1 cm.
#'@param extinction.coeff The extinction coefficient of a given protein in molar units.
#'@param mol.weight The molecular weight of a given protein. Necessary to return concentration in mg/mL.
#'
#'
#' @examples uv.vis()
#'
#' @export
uv.vis <- function(type = c("protein","dsDNA","ssDNA","ssRNA"),
                   A260, A280, dil.fac, path.length = NULL, extinction.coeff = NULL, mol.weight = NULL){

  if(is.null(path.length)){

    path = 1

    writeLines(c("----------------------------------------------------------------------------------------------",
                 "Assuming path length is 1 cm!",
                 "",
                 "If your path length differs from 1 cm, enter the correct value as 'path.length' in CENTIMETERS",
                 "-----------------------------------------------------------------------------------------------",
                 "")) }

  if(!is.null(path.length)){

    path = path.length }




  if(type == "protein"){

    writeLines(c("--------------------------------------------------------",
                 "Please enter extinction coefficient values in MOLAR units!",
                 "",
                 "Please enter protein molecular weight in GRAMS/MOLE!",
                 "--------------------------------------------------------",
                 ""))

    if(is.null(extinction.coeff)) {

      writeLines(c("Must define protein extinction coefficient in MOLAR units!",
                   "",
                   "Use the function prot.param() in this package for estimation at 280 and 205 nm.",
                   "",
                   "The web tool https://web.expasy.org/protparam/ also provides this value at 280 nm.")) }

    if(is.null(mol.weight)){

      writeLines(c("This function assumes the extinction coefficient is reported in MOLAR units!",
                   "",
                   "Molecular weight is not supplied. Only molar concentration will be reported!"))

      dat <- data.frame(A280 = A280,
                        A260 = A260,
                        A280.A260.ratio = A280/A260,
                        df = dil.fac,
                        corr.A280 = A280*dil.fac,
                        path.len = path,
                        ext.coeff = extinction.coeff)

      dat <- dat %>%
        dplyr::mutate(molar.conc = (corr.A280/(extinction.coeff*path.len))*1000000) %>%
        dplyr::select(corr.A280, molar.conc, A280.A260.ratio) %>%
        tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")


      return(dat) }

    if(!is.null(mol.weight)){

      writeLines(c("This function assumes the extinction coefficient is reported in MOLAR units!",
                   "",
                   "Concentration will be reported in both molar and mass (mg/mL) units!"))

      dat <- data.frame(A280 = A280,
                        A260 = A260,
                        A280.A260.ratio = A280/A260,
                        df = dil.fac,
                        corr.A280 = A280*dil.fac,
                        path.len = path,
                        ext.coeff = extinction.coeff)


      dat <- dat %>%
        dplyr::mutate(molar.conc = (corr.A280/(extinction.coeff*path.len))*1000000,
                      mass.conc = (molar.conc/1000000)*mol.weight) %>%
        dplyr::select(corr.A280, molar.conc, mass.conc, A280.A260.ratio) %>%
        tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

      return(dat) } }


  if(type == "dsDNA"){

    writeLines(c("Calculating concentration of dsDNA!"))

    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*50) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")


    return(dat) }

  if(type == "ssDNA"){

    writeLines(c("Calculating concentration of ssDNA!"))

    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*33) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")


    return(dat) }

  if(type == "ssRNA"){

    writeLines(c("Calculating concentration of ssRNA!"))

    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*40) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value") }


  return(dat) }
