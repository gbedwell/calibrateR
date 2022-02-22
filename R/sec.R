#' Gel filtration column calibration and analyte analysis
#'
#'
#'
#'@param void The void volume of the column, measured with something like blue dextran.
#'@param cv The column volume of the column, measured with something like 0.2% (v/v) acetone.
#'@param stds The elution volumes of the standards, entered within \code{c()}.
#'@param masses The known masses of each of the standards.
#'@param unk The elution volume of any unknown analyte whose hydrodynamic parameters are of interest. Required for any \code{type="pred"} usage.
#'
#'
#' @examples sec(masses, void, cv, stds, unk)
#'
#' @export
sec <- function(void, cv, stds, masses, unk, sample.ids) {
  require(broom)
  require(dplyr)

  std.masses <- masses
  stds <- data.frame(standards.ev = stds) %>%
    dplyr::arrange(standards.ev)
  void.vol <- void
  col.vol <- cv

  stds.df <- data.frame(mw = std.masses) %>%
    dplyr::arrange(desc(mw)) %>%
    dplyr::mutate(rh = (10^((-0.254+(0.369*log10(mw)))))/10)

  pred.ev <- data.frame(pred.ev = unk,
                        ids = sample.ids) %>%
    dplyr::mutate(ve.v0 = pred.ev - void.vol,
                  K.av = ve.v0/(max(col.vol)-void.vol)) %>%
    dplyr::arrange(K.av) %>%
    dplyr::select(ids, K.av)

  stds.df <- cbind(stds.df, stds)  %>%
    dplyr::mutate(ve.v0 = standards.ev - void.vol,
                  K.av = ve.v0/(col.vol-void.vol)) %>%
    tidyr::pivot_longer(cols = c("mw","rh"), names_to = "param", values_to = "value") %>%
    dplyr::select(-c(standards.ev, ve.v0)) %>%
    dplyr::group_by(param) %>%
    tidyr::nest(data = c(K.av, value)) %>%
    dplyr::mutate(fit = purrr::map(data, ~lm(log(value) ~ K.av, data = .)),
                  glancefit = purrr::map(fit, broom::glance),
                  tidyfit = purrr::map(fit, broom::tidy),
                  augfit = purrr::map(fit, broom::augment),
                  preds = map(fit, predict, newdata = pred.ev)) %>%
    dplyr::ungroup()

  pred.df <- stds.df %>%
    tidyr::unnest(preds) %>%
    dplyr::select(param, preds) %>%
    cbind(., pred.ev) %>%
    dplyr::select(ids, param, preds) %>%
    dplyr::mutate(preds = exp(preds)) %>%
    tidyr::pivot_wider(names_from = param, values_from = preds) %>%
    dplyr::rename(mw.da = mw,
                  rh.nm = rh)

  datlist <- list()
  datlist[[1]] <- pred.df
  datlist[[2]] <- stds.df

  return(datlist)
}
