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
#' @examples uv_vis()
#'
#' @export
uv_vis <- function(type = c("protein","dsDNA","ssDNA","ssRNA"),
                   A260, A280, dil.fac, path.length = NULL,
                   extinction.coeff = NULL, mol.weight = NULL){
  if(is.null(path.length)) {
    path = 1 }
  else {
    path = path.length }

  if(type == "protein"){
    if(is.null(extinction.coeff)) {
      stop(c("Must define protein extinction coefficient in molar units.")) }

    if(is.null(mol.weight)) {
      dat <- data.frame(A280 = A280,
                        A260 = A260,
                        A280.A260.ratio = A280/A260,
                        df = dil.fac,
                        corr.A280 = A280*dil.fac,
                        path.len = path,
                        ext.coeff = extinction.coeff)

      dat <- dat %>%
        dplyr::mutate(micromolar.conc = (corr.A280/(extinction.coeff*path.len))*1000000) %>%
        dplyr::select(corr.A280, micromolar.conc, A280.A260.ratio) %>%
        magrittr::set_colnames(c("corrected A280","concentration (µM)","280/260")) %>%
        tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

      return(dat)
    }
    else {
      dat <- data.frame(A280 = A280,
                        A260 = A260,
                        A280.A260.ratio = A280/A260,
                        df = dil.fac,
                        corr.A280 = A280*dil.fac,
                        path.len = path,
                        ext.coeff = extinction.coeff)


      dat <- dat %>%
        dplyr::mutate(micromolar.conc = (corr.A280/(extinction.coeff*path.len))*1000000,
                      mass.conc = (micromolar.conc/1000000)*mol.weight) %>%
        dplyr::select(corr.A280, micromolar.conc, mass.conc, A280.A260.ratio) %>%
        magrittr::set_colnames(c("corrected A280","concentration (µM)","concentration (mg/mL)","280/260")) %>%
        tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

      return(dat)
    }
  }

  if(type == "dsDNA"){
    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*50) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      magrittr::set_colnames(c("corrected A260","concentration (µg/mL)","260/280")) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

    return(dat)
  }

  if(type == "ssDNA"){
    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*33) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      magrittr::set_colnames(c("corrected A260","concentration (µg/mL)","260/280")) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

    return(dat)
  }

  if(type == "ssRNA"){
    dat <- data.frame(A280 = A280,
                      A260 = A260,
                      A260.A280.ratio = A260/A280,
                      df = dil.fac,
                      corr.A260 = A260*dil.fac,
                      path.len = path)

    dat <- dat %>%
      dplyr::mutate(ug.mL = corr.A260*df*path.len*40) %>%
      dplyr::select(corr.A260, ug.mL, A260.A280.ratio) %>%
      magrittr::set_colnames(c("corrected A260","concentration (µg/mL)","260/280")) %>%
      tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "value")

    return(dat) }
}
