#' Calculates standard physical prameters of a protein given its amino acid sequence.
#'
#'
#'
#'@param aa.sequence The amino acid sequence of the protein of interest in single letter code. Enter as \code{c()}.
#'@param temp The temperature of the analysis.
#'@param wavelength The wavelength being used by the refractometer.
#'
#' @examples prot_param()
#'
#' @export
prot_param <- function(aa.sequence, temp = NULL, wavelength = NULL){

  if (is.null(temp)) {
    temp = 25 }
  else { temp = temp }

  if (is.null(wavelength)) {
    wavelength = 589 }
  else { wavelength = wavelength }

  aa.sequence = stringr::str_to_upper(aa.sequence, locale = "en")

  aa.str <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  aa.dat <- data.frame(amino.acid = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
                       mono.mass = c(71.03711,156.10111,114.04293,115.02694,103.00919,129.04259,128.05858,
                                     57.02146,137.05891,113.08406,113.08406,128.09496,131.04049,147.06841,
                                     97.05276,87.03203,101.04768,186.07931,163.06333,99.06841),
                       avg.mass = c(71.0788,156.1875,114.1038,115.0886,103.1388,129.1155,128.1307,57.0519,
                                    137.1411,113.1594,113.1594,128.1741,131.1926,147.1766,97.1167,
                                    87.0782,101.1051,186.2132,163.1760,99.1326),
                       vbar = c(0.74,0.70,0.62,0.60,0.63,0.66,0.67,0.64,0.67,0.90,
                                0.90,0.82,0.75,0.77,0.76,0.63,0.70,0.74,0.71,0.86),
                       ext.h2o = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5500,1490,0),
                       ext.guhcl = c(0,0,0,0,120,0,0,0,0,0,0,0,0,0,0,0,0,5690,1280,0),
                       ext.205 = c(0,1350,400,0,690,0,400,0,5200,0,0,0,1830,8600,0,0,0,20400,6080,0),
                       dndc = c(0.167,0.206,0.192,0.197,0.206,0.183,0.186,0.175,0.219,0.179,
                                0.173,0.181,0.204,0.244,0.165,0.170,0.172,0.277,0.240,0.172))

  aa.count <- data.frame(amino.acid = aa.str,
                         aa.count = stringr::str_count(aa.sequence, aa.str))

  aa.dat <- dplyr::inner_join(aa.dat, aa.count, by="amino.acid")

  residues <- aa.dat %>%
    dplyr::summarise(value = round(sum(aa.count))) %>%
    dplyr::mutate(parameter = "Number of Residues")

  avg.mass <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*avg.mass),
                     value = round(value+18, 2)) %>%
    dplyr::mutate(parameter = "Isotopically Averaged Mass (Da)")

  mono.mass <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*mono.mass),
                     value = round(value+18, 2)) %>%
    dplyr::mutate(parameter = "Monoisotopic Mass (Da)")

  vbar.25 <- aa.dat %>%
    dplyr::mutate(parameter = "vbar at 25C (mL/g)") %>%
    dplyr::mutate(value = (sum(aa.count*(avg.mass)*vbar))/(sum(aa.count*(avg.mass)))) %>%
    dplyr::select(parameter, value) %>%
    dplyr::distinct()

  vbar.temp <- aa.dat %>%
    dplyr::mutate(parameter = "vbar at temp (mL/g)") %>%
    dplyr::mutate(vbar.25 = (sum(aa.count*(avg.mass)*vbar))/(sum(aa.count*(avg.mass))),
                  value = round(vbar.25 + (0.000425*((temp+273.15)-298.15)), 3)) %>%
    dplyr::select(parameter, value) %>%
    dplyr::distinct()

  ext.coeff <- aa.dat %>%
    dplyr::summarise(value = round(sum(aa.count*ext.h2o))) %>%
    dplyr::mutate(parameter = "Reduced Molar Ext. Coeff. at 280 nm (OD/mol*cm) in H2O")

  ext.coeff.guhcl <- aa.dat %>%
    dplyr::summarise(value = round(sum(aa.count*ext.guhcl))) %>%
    dplyr::mutate(parameter = "Unfolded Molar Ext. Coeff. at 280 nm (OD/mol*cm) in 6 M GuHCl")

  ext.205 <- aa.dat %>%
    dplyr::summarise(value = round(sum(aa.count*ext.205) + 2780*(sum(aa.count)-1))) %>%
    dplyr::mutate(parameter = "Molar Ext. Coeff. at 205 nm (OD/mol*cm) in Water")

  dndc <- aa.dat %>%
    dplyr::summarise(value.589 = (sum(aa.count*(avg.mass)*dndc))/(sum(aa.count*(avg.mass)))) %>%
    dplyr::mutate(parameter = "dn/dc (mL/g)",
                  value.578 = value.589/(0.940+(20000/(589^2))),
                  value = round(value.578*(0.940+(20000/(wavelength^2)))*(1 + (25-temp)*0.0025/30), 3)) %>%
    dplyr::select(parameter, value)

  combo <- rbind(residues, avg.mass, mono.mass, vbar.25, vbar.temp, ext.coeff, ext.205, dndc) %>%
    dplyr::select(parameter, value)

  r.sphere <- data.frame(parameter = "Radius of Sphere with Equal Mass and Density (nm)",
                         value = round((((3*(combo[2,2])*(combo[5,2]))/(4*pi*6.022E+23))^(1/3))*10^7, 2))

  r.denat <- data.frame(parameter = "Est. Denatured Radius (nm)",
                        value = round((2.21*(combo[1,2])^0.57)/10, 2))

  r.glob <- data.frame(parameter = "Est. Globular Radius (nm)",
                       value = round((4.75*(combo[1,2])^0.29)/10, 2))

  combo <- rbind(residues, avg.mass, mono.mass, vbar.temp, ext.coeff,
                 ext.coeff.guhcl, ext.205, dndc, r.sphere, r.denat, r.glob) %>%
    dplyr::select(parameter, value) %>%
    dplyr::rename(Parameter = parameter,
                  Value = value)

  return(combo)
}
