#' Calculates standard physical prameters of a protein given its amino acid sequence.
#'
#'
#'
#'@param aa.sequence The amino acid sequence of the protein of interest in single letter code. Enter as \code{c()}.
#'@param temp The temperature of the analysis. Relevant to calculation of various hydrodynamic parameters like vbar.
#'
#'
#' @examples prot.param()
#'
#' @export
prot.param <- function(aa.sequence, temp){

  aa.str <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  aa.dat <- data.frame(amino.acid = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
                       mono.mass = c(71.03711381,156.1011111,114.0429275,115.0269431,103.0091845,129.0425931,128.0585775,
                                     57.02146374,137.0589119,113.084064,113.084064,128.0949631,131.0404846,147.0684139,
                                     97.05276388,87.03202844,101.0476785,186.079313,163.0633286,99.06841395),
                       avg.mass = c(71.0779,156.18568,114.10264,115.0874,103.1429,129.11398,128.12922,57.05132,
                                    137.13928,113.15764,113.15764,128.17228,131.19606,147.17386,97.11518,
                                    87.0773,101.10388,186.2099,163.17326,99.13106),
                       vbar = c(0.74,0.70,0.62,0.60,0.63,0.67,0.66,0.64,0.67,0.90,
                                0.90,0.82,0.75,0.77,0.76,0.63,0.70,0.74,0.71,0.86),
                       ext.h2o = c(0,0,0,0,125,0,0,0,0,0,0,0,0,0,0,0,0,5500,1490,0),
                       ext.guhcl = c(0,0,0,0,120,0,0,0,0,0,0,0,0,0,0,0,0,5690,1280,0),
                       ext.205 = c(0,1350,400,0,690,0,400,0,5200,0,0,0,0,8600,0,0,0,20400,6080,0))

  aa.count <- data.frame(amino.acid = aa.str,
                         aa.count = stringr::str_count(aa.sequence, aa.str))

  aa.dat <- dplyr::inner_join(aa.dat, aa.count, by="amino.acid")

  residues <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count)) %>%
    dplyr::mutate(parameter = "number of residues")

  avg.mass <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*avg.mass),
                     value = value+18) %>%
    dplyr::mutate(parameter = "isotopically averaged mass (Da)")

  mono.mass <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*mono.mass),
                     value = value+18,
    ) %>%
    dplyr::mutate(parameter = "monoisotopic mass (Da)")

  vbar.25 <- aa.dat %>%
    dplyr::mutate(parameter = "vbar at 25C (mL/g)") %>%
    dplyr::mutate(value = (sum(aa.count*(avg.mass)*vbar))/(sum(aa.count*(avg.mass)))) %>%
    dplyr::select(parameter, value) %>%
    dplyr::distinct()

  vbar.temp <- aa.dat %>%
    dplyr::mutate(parameter = "vbar at temp (mL/g)") %>%
    dplyr::mutate(vbar.25 = (sum(aa.count*(avg.mass)*vbar))/(sum(aa.count*(avg.mass))),
                  value = vbar.25 + (0.000425*((temp+273.15)-298.15))) %>%
    dplyr::select(parameter, value) %>%
    dplyr::distinct()

  ext.coeff <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*ext.h2o)) %>%
    dplyr::mutate(parameter = "folded molar ext. coeff. at 280 nm (OD/mol*cm) in H2O")

  ext.coeff.guhcl <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*ext.guhcl)) %>%
    dplyr::mutate(parameter = "unfolded molar ext. coeff. at 280 nm (OD/mol*cm) in 6 M GuHCl")

  ext.205 <- aa.dat %>%
    dplyr::summarise(value = sum(aa.count*ext.205) + 2780*(sum(aa.count)-1)) %>%
    dplyr::mutate(parameter = "molar ext. coeff. at 205 nm (OD/mol*cm) in water")

  combo <- rbind(residues, avg.mass, mono.mass, vbar.25, vbar.temp, ext.coeff, ext.205) %>%
    dplyr::select(parameter, value)

  r.sphere <- data.frame(parameter = "radius of sphere with equal mass and density (nm)",
                         value = (((3*(combo[2,2])*(combo[5,2]))/(4*pi*6.022E+23))^(1/3))*10^7)


  r.denat <- data.frame(parameter = "denatured radius (nm)",
                        value = (2.21*(combo[1,2])^0.57)/10)


  r.glob <- data.frame(parameter = "globular radius (nm)",
                       value = (4.75*(combo[1,2])^0.29)/10)


  combo <- rbind(residues, avg.mass, mono.mass, vbar.temp,
                 ext.coeff, ext.coeff.guhcl, ext.205, r.sphere, r.denat, r.glob) %>%
    dplyr::select(parameter, value) %>%
    dplyr::mutate_if(is.numeric, round, 3)



  return(combo)

}
