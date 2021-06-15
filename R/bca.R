#' Generates a standard curve for BCA analysis and calculates the concentrations of unknown samples
#'
#'
#'
#'@param prot.conc The concentrations of the protein standards used in the assay.
#'@param blank.abs The absorbance value of the blank (buffer only).
#'@param std.abs Absorbance values of the standards at the concentrations defined in \code{prot.conc}.
#'@param abs Absorbance values of the unknown analyte(s).
#'@param dil.fac The dilution factors corresponding to each entered analyte aborbance value.
#'@param plot Defines whether or not the output should be presented graphically or not \code{"yes"} or \code{"no"}.
#'@param type Defines whether the output is column calibration or analyte analysis
#'
#'
#' @examples bca(prot.conc = c(), blank.abs, std.abs = c(), abs = NULL, plot, type)
#'
#' @export
bca <- function(prot.conc, blank.abs, std.abs, abs = NULL, dil.fac = NULL, plot = c("yes", "no"), type = c("cal", "pred")) {
  bca.df <- data.frame(conc = prot.conc,
                    abs = std.abs) %>%
    dplyr::mutate(corr.abs = abs - blank.abs)

  bca.std.fit <- lm(conc ~ corr.abs, data = bca.df)

  bca.std.r.sq <- data.frame(r.sq = summary(bca.std.fit)$r.squared)


  if(type == "cal" & plot == "no"){
    dat <- bca.df %>% dplyr::select(corr.abs)

    conc.pred <- data.frame(predict(bca.std.fit, newdata = dat, interval = "confidence")) %>%
      dplyr::mutate(conc.pred = fit,
                    conc.lwr = lwr,
                    conc.upr = upr) %>%
      dplyr::select(conc.pred, conc.lwr, conc.upr)

    bca.df <- cbind(bca.df, conc.pred, bca.std.r.sq)

    bca.df <- bca.df %>%
      mutate(across(where(is.numeric), round, digits=3))


    return(bca.df) }


  if(type == "cal" & plot == "yes"){
    dat <- bca.df %>% dplyr::select(corr.abs)

    conc.pred <- data.frame(predict(bca.std.fit, newdata = dat, interval = "confidence")) %>%
      dplyr::mutate(conc.pred = fit,
                    conc.lwr = lwr,
                    conc.upr = upr) %>%
      dplyr::select(conc.pred, conc.lwr, conc.upr)

    bca.df <- cbind(bca.df, conc.pred, bca.std.r.sq)

    plot <- ggplot(bca.df, aes(x=corr.abs, y=conc)) +
      geom_point(fill = "darkblue", shape=21, size=3, color="black") +
      geom_line(aes(x=corr.abs, y=conc.pred), color="darkblue", size=0.8) +
      geom_line(aes(x=corr.abs, y=conc.lwr), color = "black", size=0.2, linetype = "dashed") +
      geom_line(aes(x=corr.abs, y=conc.upr), color = "black", size=0.2, linetype = "dashed") +
      guides(color = FALSE,
             fill = FALSE) +
      theme_bw() +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16),
            legend.position = "top",
            legend.text = element_text(size=12),
            legend.title=element_blank(),
            legend.justification="center") +
      labs(x = "Corrected Absorbance", y = "Protein Concentration (mg/mL)")



    return(plot) }



  if(type == "pred" & (is.null(abs) | is.null(dil.fac))){

    writeLines(c("Must define both unknown analyte absorbance values and dilution factors!",
                 "",
                 "Dilution factors must be defined for every unknown analyte!")) }

  else{

    if(plot == "no") {

      bca.pred = data.frame(unk.abs = abs,
                            df = dil.fac) %>%
        dplyr::mutate(corr.abs = unk.abs - blank.abs)


      dat <- bca.pred %>% dplyr::select(corr.abs)

      conc.pred <- data.frame(predict(bca.std.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(conc.pred = fit,
                      conc.lwr = lwr,
                      conc.upr = upr) %>%
        dplyr::select(conc.pred, conc.lwr, conc.upr)

      bca.pred <- cbind(bca.pred, conc.pred) %>%
        dplyr::mutate(conc.pred.corr = conc.pred*df,
                      conc.lwr.corr = conc.lwr*df,
                      conc.upr.corr = conc.upr*df)

      bca.pred <- bca.pred %>%
        mutate(across(where(is.numeric), round, digits=3))

      return(bca.pred) }


    if(plot == "yes"){

      writeLines(c("The concentrations reported as labels on the plot are the DF CORRECTED values!",
                   "",
                   "These represent the final values that should usually be used in downstream calculations!",
                   "",
                   "For multiple dilutions of the same sample, the mean of these values should be calculated."))

      dat.std <- bca.df %>% dplyr::select(corr.abs)

      std.pred <- data.frame(predict(bca.std.fit, newdata = dat.std, interval = "confidence")) %>%
        dplyr::mutate(std.pred = fit,
                      std.lwr = lwr,
                      std.upr = upr) %>%
        dplyr::select(std.pred, std.lwr, std.upr)

      bca.df <- cbind(bca.df, std.pred, bca.std.r.sq)


      bca.pred = data.frame(unk.abs = abs,
                            df = dil.fac) %>%
        dplyr::mutate(corr.abs = unk.abs - blank.abs)


      dat <- bca.pred %>% dplyr::select(corr.abs)

      conc.pred <- data.frame(predict(bca.std.fit, newdata = dat, interval = "confidence")) %>%
        dplyr::mutate(conc.pred = fit,
                      conc.lwr = lwr,
                      conc.upr = upr) %>%
        dplyr::select(conc.pred, conc.lwr, conc.upr)

      bca.pred <- cbind(bca.pred, conc.pred) %>%
        dplyr::mutate(conc.pred.corr = conc.pred*df,
                      conc.lwr.corr = conc.lwr*df,
                      conc.upr.corr = conc.upr*df)


      plot <- ggplot(bca.df, aes(x=corr.abs, y=conc)) +
        geom_point(fill = "darkblue", shape=21, size=3, color="black") +
        geom_line(aes(x=corr.abs, y=std.pred), color="darkblue", size=0.8) +
        geom_line(aes(x=corr.abs, y=std.lwr), color = "black", size=0.2, linetype = "dashed") +
        geom_line(aes(x=corr.abs, y=std.upr), color = "black", size=0.2, linetype = "dashed") +
        geom_point(data=bca.pred, aes(x=corr.abs, y=conc.pred), shape=21, size=3, fill = "darkgreen", color="black") +
        geom_label_repel(data=bca.pred,
                        aes(x=corr.abs, y=conc.pred, label=paste(round(conc.pred.corr,2), "mg/mL")),
                        point.padding = 0.1,
                        nudge_y = 0.125,
                        color = "black") +
        guides(color = FALSE,
               fill = FALSE) +
        theme_bw() +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=16),
              legend.position = "top",
              legend.text = element_text(size=12),
              legend.title=element_blank(),
              legend.justification="center") +
        labs(x = "Corrected Absorbance", y = "Protein Concentration (mg/mL)")


    return(plot) } } }
