#' Generates a standard curve for BCA analysis and calculates the concentrations of unknown samples
#'
#'
#'
#'@param prot.conc The concentrations of the protein standards used in the assay.
#'@param blank.abs The absorbance value of the blank (buffer only).
#'@param std.abs Absorbance values of the standards at the concentrations defined in \code{prot.conc}.
#'@param abs Absorbance values of the unknown analyte(s).
#'@param dil.fac The dilution factors corresponding to each entered analyte aborbance value.
#'
#'
#' @examples bca(prot.conc = c(), blank.abs, std.abs = c(), abs = NULL)
#'
#' @export
bca <- function(prot.conc, blank.abs, std.abs, abs, sample.ids, dil.fac) {
  require(broom, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(dplyr, quietly = TRUE)

  pred.abs <- data.frame(pred.abs = abs,
                         ids = sample.ids) %>%
    dplyr::mutate(corr.abs = (pred.abs-blank.abs)*dil.fac) %>%
    dplyr::select(ids, corr.abs) %>%
    dplyr::arrange(corr.abs)

  bca.df <- data.frame(conc = prot.conc,
                       abs = std.abs) %>%
    dplyr::mutate(corr.abs = abs - blank.abs) %>%
    dplyr::select(conc, corr.abs) %>%
    tidyr::nest(data = c(conc, corr.abs)) %>%
    dplyr::mutate(fit = purrr::map(data, ~lm(conc ~ corr.abs, data = .)),
                  glancefit = purrr::map(fit, broom::glance),
                  tidyfit = purrr::map(fit, broom::tidy),
                  augfit = purrr::map(fit, broom::augment),
                  preds = map(fit, predict, newdata = pred.abs))

  pred.df <- bca.df %>%
    tidyr::unnest(preds) %>%
    cbind(., pred.abs) %>%
    dplyr::select(ids, preds) %>%
    dplyr::mutate(preds = round(preds, 3)) %>%
    dplyr::rename(`IDs` = ids,
                  `Concentration (mg/mL)` = preds)

  datlist <- list()
  datlist[[1]] <- pred.df
  datlist[[2]] <- bca.df

  return(datlist)
}
