#' Gel filtration (or size-exclusion) column calibration and analyte analysis
#'
#'
#'
#'@param void The void volume of the column, measured with something like blue dextran.
#'@param cv The column volume of the column, measured with something like 0.2% (v/v) acetone.
#'@param stds The elution volumes of the standards, entered within \code{c()}.
#'@param masses The known masses of each of the standards. If left \code{NULL}, Bio Rad standard masses will be assumed.
#'@param unk The elution volume of any unknown analyte whose hydrodynamic parameters are of interest. Required for any \code{type="pred"} usage.
#'@param plot Defines whether or not the output should be presented graphically or not \code{"yes"} or \code{"no"}.
#'@param type Defines whether the output is column calibration or analyte analysis
#'
#'
#' @examples sec(masses = c(), void, cv, stds = c(), unk = NULL, plot = "no", type = "cal")
#'
#' @export
sec <- function(void, cv, stds, masses = NULL, unk = NULL, plot = c("yes", "no"), type = c("cal", "pred")) {

  if(is.null(masses)){

    writeLines(c("No standard masses defined.",
                 "",
                 "Using Bio Rad standard masses as default.",
                 "",
                 "If these masses do not match your experiment, please define mass values.",
                 "Example: masses = c(value1, value2, etc.)"))

    stds <- data.frame(standards.ev = stds)

    void.vol = void

    col.vol = cv

    biorad.df <- data.frame(
    mw.da = c(669000,158000,44000,17000,1350)) %>%
      dplyr::mutate(rh.nm = (10^((-0.254+(0.369*log10(mw.da)))))/10)

    biorad.df <- cbind(biorad.df, stds)  %>%
      dplyr::mutate(ve.v0 = standards.ev - void.vol,
                    K.av = ve.v0/(max(col.vol)-void.vol))

    biorad.mw.fit <- lm(log(mw.da) ~ K.av, data = biorad.df)

    biorad.mw.r.sq <- data.frame(mw.r.sq = summary(biorad.mw.fit)$r.squared)

    biorad.rh.fit <- lm(log(rh.nm) ~ K.av, data= biorad.df)

    biorad.rh.r.sq <- data.frame(rh.r.sq = summary(biorad.rh.fit)$r.squared) }


  if(!is.null(masses)){

    writeLines(c("Using defined mass values instead of default values!",
                 ""))

    std.masses <- masses

    stds <- data.frame(standards.ev = stds)

    void.vol = void

    col.vol = cv

    biorad.df <- data.frame(
    mw.da = std.masses) %>%
      dplyr::mutate(rh.nm = (10^((-0.254+(0.369*log10(mw.da)))))/10)

    biorad.df <- cbind(biorad.df, stds)  %>%
      dplyr::mutate(ve.v0 = standards.ev - void.vol,
                    K.av = ve.v0/(max(col.vol)-void.vol))

    biorad.mw.fit <- lm(log(mw.da) ~ K.av, data = biorad.df)

    biorad.mw.r.sq <- data.frame(mw.r.sq = summary(biorad.mw.fit)$r.squared)

    biorad.rh.fit <- lm(log(rh.nm) ~ K.av, data= biorad.df)

    biorad.rh.r.sq <- data.frame(rh.r.sq = summary(biorad.rh.fit)$r.squared) }


  if(type == "cal" & plot == "no"){
    dat <- biorad.df %>% dplyr::select(K.av)

    mw.pred <- data.frame(predict(biorad.mw.fit, newdata = dat, interval = "confidence")) %>%
      dplyr::mutate(mw.pred = exp(fit),
                    mw.lwr = exp(lwr),
                    mw.upr = exp(upr),
                    type = "mw") %>%
      dplyr::select(mw.pred, mw.lwr, mw.upr)

    rh.pred <- data.frame(predict(biorad.rh.fit, newdata = dat, interval = "confidence")) %>%
      dplyr::mutate(rh.pred = exp(fit),
                    rh.lwr = exp(lwr),
                    rh.upr = exp(upr),
                    type = "rh") %>%
      dplyr::select(rh.pred, rh.lwr, rh.upr)

    biorad.df <- cbind(biorad.df, mw.pred, rh.pred, biorad.mw.r.sq, biorad.rh.r.sq)

    biorad.df <- biorad.df %>%
      mutate(across(where(is.numeric), round, digits=2))


    return(biorad.df) }


  if(type == "cal" & plot == "yes"){
    biorad.df <- biorad.df %>%
      tidyr::pivot_longer(cols = c(mw.da, rh.nm), names_to = "sample") %>%
      dplyr::select(sample, value, K.av) %>%
      magrittr::set_colnames(c("sample","calc.value","K.av")) %>%
      dplyr::mutate(parameter = case_when(
        sample == "mw.da"~"mw",
        sample == "rh.nm"~"rh")) %>%
     dplyr::select(-c(sample)) %>%
     dplyr::arrange(parameter) %>%
     dplyr::select(parameter, K.av, calc.value)

    dat.mw <- biorad.df %>%
      dplyr::filter(parameter == "mw") %>%
      dplyr::select(K.av)

    dat.rh <- biorad.df %>%
      dplyr::filter(parameter == "rh") %>%
      dplyr::select(K.av)

    mw.pred <- data.frame(predict(biorad.mw.fit, newdata = dat.mw, interval = "confidence")) %>%
      dplyr::mutate(pred.value = exp(fit),
                    lwr = exp(lwr),
                    upr = exp(upr)) %>%
      dplyr::select(pred.value, lwr, upr)

    rh.pred <- data.frame(predict(biorad.rh.fit, newdata = dat.rh, interval = "confidence")) %>%
      dplyr::mutate(pred.value = exp(fit),
                    lwr = exp(lwr),
                    upr = exp(upr),
                    parameter = "rh") %>%
      dplyr::select(pred.value, lwr, upr)

   pred.combine <- rbind(mw.pred, rh.pred)

   biorad.df <- cbind(biorad.df, pred.combine)

    plot <- ggplot(biorad.df, aes(x=K.av, y=log(calc.value), fill=parameter)) +
      geom_point(shape=21, size=3, color="black") +
      geom_line(aes(x=K.av, y=log(pred.value), color=parameter), size=0.8) +
      geom_line(aes(x=K.av, y=log(lwr)), color = "black", size=0.2, linetype = "dashed") +
      geom_line(aes(x=K.av, y=log(upr)), color = "black", size=0.2, linetype = "dashed") +
      scale_fill_manual(values = c("darkblue","darkred"),
                        labels = c(expression(M[W]~(Da)), expression(R[H]~(nm))),
                        guide = guide_legend(title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5)) +
     scale_color_manual(values = c("darkblue","darkred")) +
      guides(color = FALSE) +
      theme_bw() +
      theme(axis.text=element_text(size=14),
           axis.title=element_text(size=16),
           legend.position = "top",
           legend.text = element_text(size=12),
           legend.title=element_blank(),
           legend.justification="center") +
      scale_y_continuous(name = expression(log(M[W]~(Da))), limits=c(-1, 15),
                        sec.axis = sec_axis(trans=~., name=expression(log(R[H]~(nm))))) +
      scale_x_continuous(limits=c(0,1),
                         breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      xlab(expression(K[av]))



    return(plot) }


  if(type == "pred" & is.null(unk)){
    print("Must define analyte elution volumes!") }


  else{

    if(plot == "no"){

      pred.ev = data.frame(pred.ev = unk)

      biorad.pred <- pred.ev %>%
       dplyr::mutate(ve.v0 = pred.ev - void.vol,
                      K.av = ve.v0/(max(col.vol)-void.vol))

      dat <- biorad.pred %>% dplyr::select(K.av)

      mw.pred <- data.frame(predict(biorad.mw.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(mw.pred = exp(fit),
                      mw.lwr = exp(lwr),
                      mw.upr = exp(upr)) %>%
        dplyr::select(mw.pred, mw.lwr, mw.upr)

      rh.pred <- data.frame(predict(biorad.rh.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(rh.pred = exp(fit),
                      rh.lwr = exp(lwr),
                      rh.upr = exp(upr)) %>%
        dplyr::select(rh.pred, rh.lwr, rh.upr)

      biorad.pred <- cbind(biorad.pred, mw.pred, rh.pred)

      biorad.pred <- biorad.pred %>%
        mutate(across(where(is.numeric), round, digits=2))

      return(biorad.pred) }


    if(plot == "yes"){

      writeLines(c("Molecular weight and radius LABELS appearing on the plot on NOT log-transformed!",
                   "",
                   "NOTE THAT THE ACTUAL POINTS ARE STILL AT LOG-TRANSFORMED POSITIONS!",
                   "",
                   "Molecular weight LABELS are reported in kDa units.",
                   "Hydrodynamic radius LABELS are reported in nm units."))

      biorad.df <- biorad.df %>%
        tidyr::pivot_longer(cols = c(mw.da, rh.nm), names_to = "sample") %>%
        dplyr::select(sample, value, K.av) %>%
        magrittr::set_colnames(c("sample","calc.value","K.av")) %>%
        dplyr::mutate(parameter = case_when(
          sample == "mw.da"~"mw",
          sample == "rh.nm"~"rh")) %>%
        dplyr::select(-c(sample)) %>%
        dplyr::arrange(parameter) %>%
        dplyr::select(parameter, K.av, calc.value)

      dat.mw <- biorad.df %>%
        dplyr::filter(parameter == "mw") %>%
        dplyr::select(K.av)

      dat.rh <- biorad.df %>%
        dplyr::filter(parameter == "rh") %>%
        dplyr::select(K.av)

      mw.std.pred <- data.frame(predict(biorad.mw.fit, newdata = dat.mw, interval = "confidence")) %>%
        dplyr::mutate(pred.value = exp(fit),
                      lwr = exp(lwr),
                      upr = exp(upr)) %>%
        dplyr::select(pred.value, lwr, upr)

      rh.std.pred <- data.frame(predict(biorad.rh.fit, newdata = dat.rh, interval = "confidence")) %>%
        dplyr::mutate(pred.value = exp(fit),
                      lwr = exp(lwr),
                      upr = exp(upr),
                      parameter = "rh") %>%
        dplyr::select(pred.value, lwr, upr)

      pred.combine <- rbind(mw.std.pred, rh.std.pred)

      biorad.df <- cbind(biorad.df, pred.combine)


      pred.ev = data.frame(pred.ev = unk)

      biorad.pred <- pred.ev %>%
       dplyr::mutate(ve.v0 = pred.ev - void.vol,
                      K.av = ve.v0/(max(col.vol)-void.vol))

      dat <- biorad.pred %>% dplyr::select(K.av)

      mw.pred <- data.frame(predict(biorad.mw.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(mw.pred = exp(fit),
                      mw.lwr = exp(lwr),
                      mw.upr = exp(upr)) %>%
        dplyr::select(mw.pred, mw.lwr, mw.upr)

      rh.pred <- data.frame(predict(biorad.rh.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(rh.pred = exp(fit),
                      rh.lwr = exp(lwr),
                      rh.upr = exp(upr),
                      parameter = "unknown") %>%
        dplyr::select(rh.pred, rh.lwr, rh.upr, parameter)

      biorad.pred <- cbind(biorad.pred, mw.pred, rh.pred)


    plot <- ggplot(biorad.df, aes(x=K.av, y=log(calc.value), fill=parameter)) +
      geom_point(shape=21, size=3, color="black") +
      geom_line(aes(x=K.av, y=log(pred.value), color=parameter), size=0.8) +
      geom_line(aes(x=K.av, y=log(lwr)), color = "black", size=0.2, linetype = "dashed") +
      geom_line(aes(x=K.av, y=log(upr)), color = "black", size=0.2, linetype = "dashed") +
      geom_point(data=biorad.pred, aes(x=K.av, y=log(mw.pred)), shape=21, size=3, color="black") +
      geom_point(data=biorad.pred, aes(x=K.av, y=log(rh.pred)), shape=21, size=3, color="black") +
      geom_label_repel(data=biorad.pred,
                        aes(x=K.av, y=log(mw.pred), label=round(mw.pred/1000,2)),
                        point.padding = 0.2,
                        nudge_y = 0.125,
                        color = "black",
                       fill="white") +
      geom_label_repel(data=biorad.pred,
                        aes(x=K.av, y=log(rh.pred), label=round(rh.pred,2)),
                        point.padding = 0.2,
                        nudge_y = 0.125,
                        color = "black",
                       fill="white") +
      scale_fill_manual(values = c("darkblue","darkred","darkgreen"),
                        labels = c(expression(M[W]~(Da)), expression(R[H]~(nm)), expression(Unknown)),
                        guide = guide_legend(title.position = "top",
                                             nrow = 1,
                                             title.hjust=0.5)) +
      scale_color_manual(values = c("darkblue","darkred","darkgreen")) +
      guides(color = FALSE) +
      theme_bw() +
      theme(axis.text=element_text(size=14),
           axis.title=element_text(size=16),
           legend.position = "top",
           legend.text = element_text(size=12),
           legend.title=element_blank(),
           legend.justification="center") +
      scale_y_continuous(name = expression(log(M[W]~(Da))), limits=c(-1, 15),
                        sec.axis = sec_axis(trans=~., name=expression(log(R[H]~(nm))))) +
      scale_x_continuous(limits=c(0,1),
                         breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      xlab(expression(K[av]))


    return(plot) } } }
