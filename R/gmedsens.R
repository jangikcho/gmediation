#' Probabilistic Sensitivity Analysis for Generalized Causal Mediation
#'
#' Extends the mediation/path effect framework of \code{\link{gmediate}} with
#' a probabilistic sensitivity analysis. Path effects that involve
#' non-identifiable components (paths whose identification requires an
#' assumption about the correlation between residuals of two mediators, or
#' between a mediator and the outcome, that cannot be estimated from the
#' observed data) are re-evaluated across a range of assumed correlations, so
#' that the robustness of the substantive conclusions to this un-testable
#' assumption can be assessed.
#'
#' As in \code{\link{gmediate}}, the user supplies already-fit models for the
#' outcome and for each first- and second-stage mediator; these are re-fit
#' internally on the complete-case data, and bootstrap resampling is used to
#' construct confidence intervals. \code{gmedsens} additionally varies the
#' assumed correlation(s) between residuals for non-identifiable paths over
#' the grid of values given in \code{xcorrv}, and summarizes how the path
#' estimates and their confidence intervals change across that grid, both in
#' tabular form and in sensitivity plots (rendered using \pkg{ggplot2} and
#' \pkg{gridExtra}).
#'
#' As with \code{\link{gmediate}}, the default reference group, \code{ref =
#' "entire sample"}, uses the whole complete-case sample; a different
#' reference group may instead be specified as a vector giving the reference
#' group variable name and value(s).
#'
#' \code{gmedsens} is called primarily for its side effects: it prints a
#' summary of the analysis and sensitivity tables/plots to the console /
#' graphics device, when \code{plot_out[1] = "true"} writes all sensitivity
#' plots to a single PDF in the directory given in \code{plot_out[2]}, and,
#' when \code{table_out[1] = "true"}, writes the resulting path-specific,
#' first-stage, and second-stage sensitivity tables as CSV files to the
#' directory given in \code{table_out[2]}.
#'
#' @param data a data frame containing all variables used in the outcome and
#'   mediator models.
#' @param model.m1 a list of the model fit(s) for the first-stage mediator(s).
#' @param model.m2 a list of the model fit(s) for the second-stage
#'   mediator(s).
#' @param model.y an object with the model fit for the response (outcome)
#'   variable.
#' @param expos a vector with the name of the exposure variable and the
#'   exposure value, e.g. \code{c("exposure variable name", exposure value)}.
#' @param ref reference group specification. Either the default
#'   \code{"entire sample"}, meaning the whole complete-case sample is used
#'   as the reference group, or a vector with the name of the reference group
#'   variable and the reference group value(s), e.g. \code{c("reference group
#'   variable name", c(reference group value))}.
#' @param refmult numeric multiplier used to increase the number of Monte
#'   Carlo draws taken from the reference group, to reduce Monte Carlo error.
#' @param seed seed used to generate bootstrap samples.
#' @param boots number of bootstrap samples to generate.
#' @param cluster the cluster variable used for cluster bootstrap resampling,
#'   if applicable. Default is \code{NULL} (individual-level resampling).
#' @param decomp integer (1 or 2) specifying which of two available
#'   decompositions of the total exposure effect is used to define path
#'   effects.
#' @param CIPCT confidence level (as a percentage) used for the bootstrap
#'   confidence intervals.
#' @param xcorrv a vector of sensitivity parameter values (assumed
#'   correlations between residuals) evaluated for non-identifiable paths.
#' @param scale character string specifying the scale on which path effects
#'   are summarized (e.g. \code{"MeanDiff"}).
#' @param deltsi numeric tolerance used when checking whether a computed
#'   correlation/covariance falls within a valid range.
#' @param tol numeric convergence tolerance used in the root-finding /
#'   probability-matching steps of the sensitivity calculations.
#' @param tolx numeric convergence tolerance used for the argument (x) in the
#'   root-finding steps of the sensitivity calculations.
#' @param assumvars character vector naming the mediator/outcome variables
#'   for which non-identifiable paths should be evaluated across
#'   \code{xcorrv}; \code{"ALL"} (the default) evaluates all applicable
#'   variables.
#' @param g numeric factor applied to scale the assumed between-mediator (or
#'   mediator-outcome) correlation grid.
#' @param n_rows number of rows of panels used when laying out sensitivity
#'   plots on a page.
#' @param n_cols number of columns of panels used when laying out sensitivity
#'   plots on a page.
#' @param ndig number of digits used when rounding values displayed in the
#'   sensitivity tables.
#' @param rows_per_page maximum number of table rows to display per printed
#'   page/table when results span many rows.
#' @param plot_out a length-2 character vector controlling whether the
#'   sensitivity plots are written to disk. The first element is
#'   \code{"true"} or \code{"false"} (the default), and the second is the
#'   directory in which to write the PDF \code{gmedsens_all_plots.pdf}
#'   (created if it does not exist) when the first element is \code{"true"}.
#'   The PDF contains the plots for every pair of assumption variables
#'   evaluated. If the Cairo graphics device is unavailable, a standard PDF
#'   device is used instead (in which case arrow glyphs in plot titles may
#'   render as dots).
#' @param table_out a length-2 character vector controlling whether the
#'   sensitivity summary tables are written to disk. The first element is
#'   \code{"true"} or \code{"false"} (the default), and the second is the
#'   directory (created if it does not exist) in which to write the CSV files
#'   \code{gmedsens_pathspec_out.csv}, \code{gmedsens_firststage_out.csv} and
#'   \code{gmedsens_secondstage_out.csv}. When several pairs of assumption
#'   variables are evaluated (\code{assumvars = "ALL"}), the files contain the
#'   results for every pair, identified by the \code{Rho.Between} column.
#' @param maxcnt maximum number of iterations allowed in the internal
#'   root-finding routines used for the sensitivity calculations.
#'
#' @return Called primarily for its side effects (see Details); returns
#'   \code{NULL} invisibly.
#'
#' @seealso \code{\link{gmediate}} for the underlying causal mediation and
#'   path analysis this function builds on.
#'
#' @examples
#' \donttest{
#' model.y2 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + visit +
#'                    seal + ohi, family = binomial, data = dental)
#' model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
#'                   family = binomial, data = dental)
#' model.m12 <- glm(visit ~ ses + grpvlb + grpbpd + race + SEX,
#'                   family = binomial, data = dental)
#' model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
#'                   family = binomial, data = dental)
#' model.m22 <- glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
#'                   family = gaussian, data = dental)
#'
#' gmedsens(data = dental, model.m1 = list(model.m11, model.m12),
#'          model.m2 = list(model.m21, model.m22), model.y = model.y2,
#'          expos = c("ses", 1), seed = 443095, boots = 100,
#'          cluster = "momid", decomp = 2, CIPCT = 95,
#'          xcorrv = c(0, 0.5, 0.9), scale = "MeanDiff", assumvars = c("dmftd",
#'          "ohi"), g = 1, ndig = 4, n_rows = 1, n_cols = 1,
#'          rows_per_page = 18, maxcnt = 100,
#'          plot_out = c("true", tempdir()), table_out = c("true", tempdir()))
#' }
#'
#' @importFrom stats na.omit predict quantile ecdf getCall pnorm pbinom ppois
#'   pnbinom rbinom rnbinom rnorm rpois runif dnorm qnorm complete.cases
#' @importFrom MASS glm.nb rnegbin
#' @importFrom plyr count
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_ribbon labs theme_bw theme
#'   scale_x_continuous element_text
#' @importFrom grid textGrob gpar unit grid.newpage grid.draw
#' @importFrom gridExtra tableGrob grid.arrange arrangeGrob
#' @importFrom glue glue
#' @importFrom utils write.csv tail
#' @importFrom VGAM probitlink
#' @importFrom dplyr select %>%
#' @export
gmedsens = function(data, model.m1, model.m2, model.y, expos, ref = "entire sample",
                    refmult = 10, seed = 443095, boots = 499, cluster = NULL,
                    decomp = 2, CIPCT = 95, xcorrv = c(-0.9, 0, 0.5, 0.9), 
                    scale = "MeanDiff", deltsi = 0.05, tol = 0.0001, tolx = 0.0001, 
                    assumvars = "ALL", g = 1, n_rows = 1, n_cols = 1, ndig = 4,
                    rows_per_page = 16, plot_out = c("false",""), table_out = c("false",""), maxcnt = 100){

  #############################################
  ## Count number of mediators in each stage ##
  #############################################
  if(class(model.m1)[1] == "list"){
    if(is.null(names(model.m1))){
      numberofm1 = length(model.m1)
    } else {
      cc = count(names(model.m1) == "y")
      numberofm1 = cc$freq[cc$x == TRUE]
    }
  } else {
    numberofm1 = 1
  }
  
  if(class(model.m2)[1] == "list"){
    if(is.null(names(model.m2))){
      numberofm2 = length(model.m2)
    } else {
      cc = count(names(model.m2) == "y")
      numberofm2 = cc$freq[cc$x == TRUE]
    }
  } else {
    numberofm2 = 1
  }
  
  #############################
  ## Define Warning Messages ##
  #############################
  if("th.warn" %in% names(model.y)){
    stop("please check your model convergence")
  }
  
  if(class(model.m1)[1] == "list"){
    if(is.null(names(model.m1))){
      for (i in numberofm1){
        if ("th.warn" %in% names(model.m1[[i]])){
          stop("please check your model convergence")
          
        }
      }
    } else{
      if ("th.warn" %in% names(model.m1)){
        stop("please check your model convergence")
      }
    }
  } else {
    if ("th.warn" %in% names(model.m1)){
      stop("please check your model convergence")
    }
  }
  
  if(class(model.m2)[1] == "list"){
    if(is.null(names(model.m2))){
      for (i in numberofm2){
        if ("th.warn" %in% names(model.m2[[i]])){
          stop("please check your model convergence")
        }
      }
    } else{
      if ("th.warn" %in% names(model.m2)){
        stop("please check your model convergence")
      }
    }
  } else {
    if ("th.warn" %in% names(model.m2)){
      stop("please check your model convergence")
    }
  }
  
  # Ensure g is within range (0, 1]; default to 1 if out of range
  if (g <= 0 | g > 1) {
  warning("Invalid g value. g must be in the range (0, 1]. Defaulting to g = 1.")
  g <- 1
}
  
  ##############
  ## Set seed ##
  ##############
  set.seed(seed)
  
  #########################################
  ## Select complete cases from the data ##
  #########################################
  gdata = na.omit(data)
  
  ######################################################
  ## Identify response and covariates from the models ##
  ######################################################
  # model y
  yvar = names(model.y$model)
  y = yvar[1]
  y.allcov = yvar[-1]
  
  # model m2
  m2 = vector("list", numberofm2)
  m2.eachcov = vector("list", numberofm2)
  
  for (i in 1:numberofm2){
    m2[[i]] = names(model.m2[[i]]$model)[1]
    m2.eachcov[[i]] = names(model.m2[[i]]$model)[-1]
  }
  
  # model m1
  m1 = vector("list", numberofm1)
  m1.eachcov = vector("list", numberofm1)
  
  for (i in 1:numberofm1){
    m1[[i]] = names(model.m1[[i]]$model)[1]
    m1.eachcov[[i]] = names(model.m1[[i]]$model)[-1]
  }
  
  # all covariates in m2 models
  m2.tempcov = c()
  for (i in 1:numberofm2){
    m2.tempcov = c(m2.tempcov, m2.eachcov[[i]])
  }
  m2.allcov = unique(m2.tempcov)
  
  # all covariates in m1 models
  m1.tempcov = c()
  for (i in 1:numberofm1){
    m1.tempcov = c(m1.tempcov, m1.eachcov[[i]])
  }
  m1.allcov = unique(m1.tempcov)
  
  ##-- choose non overlapping covariates with m1 and m2 from y model
  overlap.y = which(y.allcov %in% c(m2,m1))
  if(length(overlap.y) == 0){
    y.cov = y.allcov
  } else {
    y.cov = y.allcov[-overlap.y]
  }
  
  ##-- choose non overlapping covariates with m1 from m2 model
  overlap.m2 = which(m2.allcov %in% c(m1))
  if(length(overlap.m2) == 0){
    m2.cov = m2.allcov
  } else {
    m2.cov = m2.allcov[-overlap.m2]
  }
  
  m1.cov = m1.allcov
  
  covariates = c(unique(c(y.cov, m2.cov, m1.cov)))
  
  m1aug = c("   " , m1)  
  m2aug = c("   " , m2)  
  
  #################################################################################################################
  ## y : y outcome in vector format                                                                              ##
  ## m2 : second stage mediators in list format                                                                  ##
  ## m1 : first stage mediators in list format                                                                   ##
  ## y.allcov : independent variables in y model                                                                 ##
  ## m2.eachcov : indepedent variables in each models in list format                                             ##
  ## m1.eachcov : indepednet variables in each m1 models in list format                                          ##
  ## m2.allcov : union of independent variables from m2 models in vector format                                  ##
  ## m1.allcov : union of independent variables from m1 models in vector format                                  ##
  ## y.cov : independent variables in y model excluding all second and first stage mediators in vector format    ##
  ## m2.cov : independent variables in m2 models excluding all first stage mediators in vector format            ##
  ## m1.cov : independent variables in m1 models                                                                 ##
  ## covariates : union of all independent variables excluding first and second stage mediators in vector format ##
  #################################################################################################################
  
  #############################################
  ## Check for required model specifications ##
  #############################################
  ##- define initial functions and values
  # negate function
  `%ni%` = Negate(`%in%`)
  
  # number of possible pathways
  npath = (numberofm1 + 1) * (numberofm2 + 1)
  
  # number of specified x-world correlations
  nxcorr = length(xcorrv)
  
  # NEW - NUmber of sens param (rhosi) values
  nrhosi =  2 / deltsi       # use delta as input (needs to result in whole number)
  nrcnt = nrhosi - 1
   
  #  using 2 above since go from -1 to 1
  ##- check model specifications
  # need at least 1 1st stage and 1 2nd stage med
  if (numberofm1 * numberofm2 == 0) {
    stop("Need at least 1 mediator for each stage")
  }
  
  # check second stage mediator model specifications
  for (k in 1:numberofm2) {
    # check if M2s are not exposure
    if (expos[1] == m2[[k]]) {
      stop("Exposure cannot be mediator")
    }
    # check if Y is not M2s
    if (y == m2[[k]]) {
      stop("Final outcome cannot be mediator")
    }
    # check if at least 1 1st stage mediator in each m2 model
    check2 = FALSE
    for (j in 1:numberofm1) {
      if (m1[[j]] %in% m2.eachcov[[k]]) {
        check2 = TRUE
        if (check2) break
      }
    }
    if (!check2) {
      stop("At least one M1 must be in each M2 model")
    }
    
    # check if each M2 in Y model
    if (m2[[k]] %ni% y.allcov) {
      stop("Each M2 must be in Y model")
    }
  }
  
  # check first stage mediator model specifications
  for (j in 1:numberofm1) {
    # check if exposure is in every M1 model
    if (expos[1] %ni% m1.eachcov[[j]]) {
      stop("Exposure must be in each M1 model")
    }
    # check if M1s are not exposure
    if (expos[1] == m1[[j]]) {
      stop("Exposure cannot be mediator")
    }
    # check if Y is not M1s
    if (y == m1[[j]]) {
      stop("Final outcome cannot be mediator")
    }
  }
  
  # check that each M1 is in some M2 or in the Y model
  for (j in 1:numberofm1) {
    check12 = FALSE
    for (k in 1:numberofm2) {
      if (m1[[j]] %in% m2.eachcov[[k]]) {
        check12 = TRUE
        break
      }
    }
    if (m1[[j]] %in% y.cov) {
      check12 = TRUE
    }
    if (check12 == FALSE) {
      stop("Each M1 must be a covariate in the Y or some M2 model")
    }
  }
  
  #################################
  ##- Collect Model Information -##
  #################################
  # distribution family
  family.y = list()
  family.m2 = list()
  family.m1 = list()
  
  # call for model
  call.y = list()
  call.m2 = list()
  call.m1 = list()
  
  # obs predictors for models 
  ycovars = list()
  m2covars = list()
  
  # variables in each model
  names.y = list()
  names.m2 = list()
  names.m1 = list()
  
  # Initialize "variables" matrix using variables used in Y model
  variables = matrix(names(model.y$model))
  
  # model information from y model
  family.y = model.y$family$family
  isNegBin = grep("Negative Binomial", substr(family.y,1,17))
  if (length(isNegBin) == 1) {family.y = "negbin"}
  
  call.y = getCall(model.y)
  
  # model information from 2nd stage mediator's models
  for (k in 1:numberofm2) {
    # model can be from glm or glm.nb
    # rename Negative Binomial from glm.nb to negbin
    family.m2[[k]] = model.m2[[k]]$family$family
    isNegBin = grep("Negative Binomial", substr(family.m2[[k]],1,17))
    if (length(isNegBin) == 1) {family.m2[[k]] = "negbin"}
    
    # call of model for each m2
    call.m2[[k]] = getCall(model.m2[[k]])
    
    names.m2[[k]] = matrix(names(model.m2[[k]]$model))
    variables = rbind(variables ,names.m2[[k]])
  }
  
  # model information from 1st stage mediator's models
  for (j in 1:numberofm1) {
    # model can be from glm or glm.nb
    # rename Negative Binomial from glm.nb to negbin
    family.m1[[j]] = model.m1[[j]]$family$family
    isNegBin = grep("Negative Binomial", substr(family.m1[[j]],1,17))
    if (length(isNegBin == 1)) {family.m1[[j]] = "negbin"}
    
    # call of model for each m1
    call.m1[[j]] = getCall(model.m1[[j]])
    
    names.m1[[j]] = matrix(names(model.m1[[j]]$model))
    variables = rbind(variables, names.m1[[j]])
  }
  
  if (length(cluster) == 1) {
    variables = rbind(variables, matrix(cluster))
  }
  
  ###################################################################################
  ##-- family.y : list of distribution family from y model                       --##
  ##-- family.m2 : list of distribution families from 2nd stage mediators models --##
  ##-- family.m1 : list of distribution families from 1st stage mediators models --##
  ##-- call.m2 : model information from 2nd stage mediators models               --##
  ##-- call.m1 : model information from 1st stage mediators models               --##
  ##-- variables : variables used in the entire models                           --##
  ###################################################################################
  
  ###############################################
  ##-- data.cr : dataset for reference group --##
  ###############################################
  # data set for reference group
  pred.dat.y = list()
  pred.dat.m2 = list()
  
  # m2s for opposite ds
  pred.dat.m2o = list() 
  
  pred.dat.m1.c = list()
  pred.dat.m1.e = list()
  
  ##-- Define reference group
  if (identical(ref, "entire sample")){
    data.cr0 = gdata
  } else {
    data.cr0 = c()
    for (i in 2:length(ref)){
      data.cr0 = rbind(data.cr0, gdata[gdata[,ref[1]] == ref[i],])
    }
  } 
    data.cr = c()
    if(refmult == 1){
      data.cr = data.cr0
    } else {
      for (i in 1:refmult){
        data.cr = rbind(data.cr0, data.cr)
      }
    }

  ##-- Define dataset for the reference group in Y
  pred.dat.y = data.cr[,y.allcov]
  
  ##-- Define dataset for the reference group in the 2nd stage mediators
  for (k in 1:numberofm2) {
    pred.dat.m2[[k]] = data.cr[,m2.eachcov[[k]]]
    # m2 predictors for opposite values for ds
    pred.dat.m2o[[k]] = data.cr[,m2.eachcov[[k]]]
  }
  
  ##-- Define dataset for the reference group in the 1st stage mediators and assume all are either control or exposed
  lexp = levels(as.factor(gdata[,expos[1]]))
  exposed = as.numeric(expos[2])
  control = as.numeric(lexp[lexp != expos[2]])
  
  for (j in 1:numberofm1) {
    pred.dat.m1.e[[j]] = pred.dat.m1.c[[j]] = data.frame(data.cr[,m1.eachcov[[j]]])
    colnames(pred.dat.m1.e[[j]]) = c(m1.eachcov[[j]])
    colnames(pred.dat.m1.c[[j]]) = c(m1.eachcov[[j]])
    pred.dat.m1.e[[j]][,expos[1]] = exposed
    pred.dat.m1.c[[j]][,expos[1]] = control
  }
  
  #############################################################################################################
  ##-- pred.dat.y : dataset for Y in the reference group                                                   --##
  ##-- pred.dat.m2 : dataset for 2nd stage mediators in the reference group                                --##
  ##-- pred.dat.m1.e : dataset for 1st stage mediators in the reference group assuming everyone is exposed --##
  ##-- pred.dat.m1.c : dataset for 1st stage mediators in the reference group assuming everyone is control --##
  #############################################################################################################
  
  ####################################################
  ##-- Initialize matrix to define decompositions --##
  ####################################################
  indexDM = matrix(NA, 2, npath)
  cnt = 1
  
  ##-- Get decomposition order for decomp = 1 (reverse for decomp = 2)
  for (j in 1:(numberofm1 + 1)) {
    for (k in 1:(numberofm2 + 1)) {
      indexDM[,cnt] = c(j,k)
      cnt = cnt + 1
    }
  }
  
  if (decomp == 2) {
    indexDM = indexDM[,npath:1]
  }
  
  ##################################################################
  ### Loop of assumptions - i1, i2 pairs (do sens anal for each) ###
  ##################################################################
  
  if (all(assumvars != "ALL")) {

      if (assumvars[1] == y) {
        i1 = 1
        for (i in 1:numberofm1){
          if (assumvars[2] == m1[[i]]) {i2 = i}
        }
        for (j in 1:numberofm2){
          if (assumvars[2] == m2[[j]]) {i2 = numberofm1 +j}
        }
      }
      else { 
        for (j in 1:numberofm2){
          if (assumvars[1] == m2[[j]]) {i1 = j+1}
        }
        for (i in 1:numberofm1){
          if (assumvars[2] == m1[[i]]) {i2 = i}
        }
      }
    avar1 = c(i1:i1)
    avar2 = c(i2:i2)
  }
  
  ##################################################################
  ### Combined outputs shared by ALL assumption pairs (i1, i2)   ###
  ##################################################################
  
  all_df  <- list()   # path-specific tables, one element per pair
  all_df1 <- list()   # first-stage mediator tables, one element per pair
  all_df2 <- list()   # second-stage mediator tables, one element per pair
  pdf_dev <- NULL     # PDF device; only opened when plot_out[1] == "true"
  out_plot <- NULL
  page_number <- 1
  
  if (plot_out[1] == "true") {
    if (nchar(plot_out[2]) == 0) {
      cat("WARNING: Please Specify a Folder in plot_out[2] to Save the Plots PDF\n")
    } else {
      
      # Create the output folder if it does not exist yet
      if (!dir.exists(plot_out[2])) {
        dir.create(plot_out[2], recursive = TRUE, showWarnings = FALSE)
      }
      
      out_plot <- file.path(plot_out[2], "gmedsens_all_plots.pdf")
      
      # Prefer cairo_pdf() (renders the arrows and rho correctly). On some Macs it
      # only WARNS ("failed to load cairo DLL") and opens no device, so treat a
      # warning, an error, or "no new device" all as failure and fall back to the
      # base pdf() device so a file is still produced.
      devs_before <- grDevices::dev.list()
      dev_open <- tryCatch({
        grDevices::cairo_pdf(out_plot, width = 8, height = 6, onefile = TRUE)
        length(grDevices::dev.list()) > length(devs_before)
      }, error = function(e) {
        cat("cairo_pdf() failed:", conditionMessage(e), "\n")
        FALSE
      }, warning = function(w) {
        cat("cairo_pdf() unavailable:", conditionMessage(w), "\n")
        FALSE
      })
      if (!dev_open) {
        cat("Using base pdf() instead of cairo_pdf(); arrow symbols in path labels will show as dots.\n")
        grDevices::pdf(out_plot, width = 8, height = 6, onefile = TRUE)
      }
      pdf_dev <- grDevices::dev.cur()
      
      # Close the PDF device if the function stops early (error or interrupt)
      on.exit({
        if (!is.null(pdf_dev) && pdf_dev %in% grDevices::dev.list()) grDevices::dev.off(pdf_dev)
      }, add = TRUE)
    }
  }
  
if (all(assumvars == "ALL")) {avar1 = c(1:(numberofm2 + 1))}
  
  for (i1 in avar1)  {
    if (all(assumvars == "ALL")) {
      if (i1==1) {
        i2end = numberofm1 + numberofm2}
      if (i1 > 1) {
        i2end = numberofm1}
      avar2 = c(1:i2end)
    }

  for (i2 in avar2)  {
    
  ####################################
  ##-- Generate Bootstrap Samples --##
  ####################################
  bootsamples = vector("list", boots + 1)
  
  if (is.null(cluster)){
    nclus = nrow(gdata)
    gdata$numsamp = c(1:nclus)
    
    for (i in 1:boots){
      samp = sample(gdata$numsamp, size = nclus, replace = TRUE)
      bootsamples[[i]] = gdata[samp,]
    }
    bootsamples[[boots+1]] = gdata
  }
  
  if (!is.null(cluster)){
    nclus = length(unique(gdata[,cluster]))
    
    for (i in 1:boots){
      samp = sample(gdata[,cluster], size = nclus, replace = TRUE)
      bootsamples[[i]] = data.frame()
      for (j in 1:length(samp)){
        bootsamples[[i]] = rbind(bootsamples[[i]],gdata[gdata[,cluster] == samp[j],])
      }
    }
    bootsamples[[boots+1]] = gdata
  }
  
  ##########################
  ##-- Initialize lists --##
  ##########################
  modelb.y = c()
  modelb.m2 = list()
  modelb.m1 = list()
  
  mu.m2 = list()
  mu.m2o = list()  
  
  mu.m2c = list()

  mu.m1 = list()
  mu.m12 = list()
  mu.m1.e = list()
  mu.m1.c = list()
  mu.m1.d00 = list()  
  mu.m1.dj0 = list()  
  mu.m1.d0k = list()  
  
  mu.m1.djk = list()  
  
  pred.m2 = list()
  pred.m12 = list()
  pred.m1.e = list()
  pred.m1.c = list()
  pred.m1 = list()
  PredM12V = vector("list", numberofm1)
  
  pred.m1c = list()  
  disp.m2 = list()    # list with dispersion par values 2nd stage
  disp.m1 = list()    # list with dispersion par values 1st stage
  
  # std devs for corrections
  # covers both normal and binomial cases
  stdv.m2 = list()    # list with dispersion par values 2nd stage
  stdv.m2c = list()   
  stdv.m1 = list()    # list with dispersion par values 1st stage
  
  xcrep = list()
  xcrepd = list()
  nonnullm1 = list() 
  IDENT = list()    
  
  #initialize nonnullm1
  
  for (i in 1:numberofm1) {
    nonnullm1[[i]] = matrix(0,1,numberofm2+1)
  }
  
  EYDM = list()      # exp pot outcome matrices for boots, rhos
  EDM = list()      # expected path effects

  EDM_SC = list()    # expected path effects for different scales
  EDMR = list()     # expected path effects ratios 

  PV = vector("list", npath)       # P values
  PV_scale = vector("list", npath)
  
  EDM1 = vector("list", numberofm1)    
  for (t in 1:numberofm1){
    EDM1[[t]] = matrix(0, boots+1, nrhosi-1)
  }
  
  EDM1_SC = vector("list", numberofm1)   
  for (t in 1:numberofm1){
    EDM1_SC[[t]] = matrix(0, boots+1, nrhosi-1)
  }
  
  EDM2 = vector("list", numberofm2)  
  for (t in 1:numberofm2){
    EDM2[[t]] = array(0, c(boots+1, nrhosi-1, nxcorr))
  }
  
  EDM2_SC = vector("list", numberofm2)  
  for (t in 1:numberofm2){
    EDM2_SC[[t]] = array(0, c(boots+1, nrhosi-1, nxcorr))    
  }
  
  EDT = matrix(0, boots+1, nrhosi-1)
  EDT_SC = matrix(0, boots+1, nrhosi-1)  
  
  xcorr = c()
  
  pathpr = list()   
  
  ########################
  ##-- Bootstrap loop --##
  ########################
  
  #initialize for each EY
  for (p in 1:(npath+1)) {
     EYDM[[p]] = array(0,c(boots+1,nrhosi-1,nxcorr))
  }
 
  for (c in 1 : npath){
      EDM[[c]] = array(0,c(boots+1,nrhosi-1,nxcorr))  
  }
  
  for (c in 1 : npath){
     EDM_SC[[c]] = array(0,c(boots+1,nrhosi-1,nxcorr)) 
  }
  
  for (c in 1 : npath){
    EDMR[[c]] = matrix(rep(0,((boots+1)*(nrhosi-1))),boots+1,nrhosi-1)
  }
  
  EDMM1 = array(0,c(npath,nrcnt,boots+1))
  EDMM1_SC = array(0,c(npath,nrcnt,boots+1))  
  
  EDMM2 = array(0,c(npath,nrcnt,boots+1,nxcorr))
  EDMM2_SC = array(0,c(npath,nrcnt,boots+1,nxcorr))   
  
  for (b in 1:(boots + 1)){
    
    nullv = matrix(0,1,npath)
    DISC = matrix(0,1,npath+1)
    
    ## Fit models for each sample (boot and orig)
    for (k in 1:numberofm2) {
      call.m2[[k]]$data = bootsamples[[b]][, c(m2[[k]], m2.eachcov[[k]])]
      m2covars[[k]] = bootsamples[[b]][, c(m2.eachcov[[k]])]
      
      if (family.m2[[k]] == "negbin"){
        modelb.m2[[k]] = glm.nb(call.m2[[k]]$formula, call.m2[[k]]$data)
      } else {
        modelb.m2[[k]] = eval.parent(call.m2[[k]])
        
        # Compute estim std devs - crude marginal prop for now
        if (family.m2[[k]] == "binomial"){
          pm2 = mean(bootsamples[[b]][, c(m2[[k]])])
          stdv.m2[[k]] = sqrt(pm2 * (1-pm2))
        }
      }
      
      if (family.m2[[k]] == "gaussian") {
        disp.m2[[k]] = summary(modelb.m2[[k]])$disp
        #stdv.m2[[k]] = disp.m2[[k]]
        stdv.m2[[k]] = sqrt(disp.m2[[k]])  # JA 7-28-25
      }
      
      if (family.m2[[k]] == "negbin") {
        disp.m2[[k]] = summary(modelb.m2[[k]])$theta
        stdv.m2[[k]] = disp.m2[[k]]
      }
      
      # Get means for indivs in ref grp
      mu.m2[[k]] = predict(modelb.m2[[k]], type = "response", newdata = pred.dat.m2[[k]])
      
      if (family.m2[[k]] == "poisson") {
        stdv.m2[[k]] = sqrt(mu.m2[[k]])
      }
    }
    
    for (j in 1:numberofm1){
      # Get data (all model vars) for bootstrap sample for refit
      
      # should be the individual indices acc to cluster boot draw
      call.m1[[j]]$data = bootsamples[[b]][, c(m1[[j]], m1.eachcov[[j]])]
      
      # Refit model using current sample (boot or orig)
      if (family.m1[[j]] == "negbin"){
        modelb.m1[[j]] = glm.nb(call.m1[[j]]$formula, call.m1[[j]]$data)
      } else {
        modelb.m1[[j]] = eval.parent(call.m1[[j]])
        
        if (family.m1[[j]] == "binomial"){
          pm1 = mean(bootsamples[[b]][, c(m1[[j]])])
          stdv.m1[[j]] = sqrt(pm1 * (1-pm1))
        }
      }
      
      if (family.m1[[j]] == "gaussian" | family.m1[[j]] == "negbin") {
        disp.m1[[j]] = summary(modelb.m1[[j]])$disp
        stdv.m1[[j]] = sqrt(disp.m1[[j]])  
      }
      
      # Get means for indivs in ref grp under each expos status
      mu.m1.e[[j]] = predict(modelb.m1[[j]], type = "response", newdata = pred.dat.m1.e[[j]])
      mu.m1.c[[j]] = predict(modelb.m1[[j]], type = "response", newdata = pred.dat.m1.c[[j]])
    }
    
    call.y$data = bootsamples[[b]][, c(y, y.allcov)]
    ycovars = bootsamples[[b]][ , c(y.allcov)]     
    
    if (family.y == "negbin"){
      modelb.y = glm.nb(call.y$formula, call.y$data)
      disp.y = summary(modelb.y)$theta
    } else {
      modelb.y = eval.parent(call.y)    # refit Y model

      if (family.y == "gaussian"){
         disp.y = summary(modelb.y)$disp
         stdv.y = sqrt(disp.y)  
      }      

      if (family.y == "binomial"){
         py = mean(bootsamples[[b]][, c(y)])
         stdv.y = sqrt(py * (1-py))
         }
    }
    
    ##-- End of fitting models
    
    ###################################################################
    ###    Med Formula: Draw values for each med, each path eff     ###
    ###################################################################
    
    DM = matrix(1, numberofm1+1, numberofm2+1)
    
    first = matrix(1, numberofm1, nxcorr)
    
    nref = dim(data.cr)[1]
    
    # get pred values for each m1, each expos status
    
    for (j in 1:numberofm1) {
      if (family.m1[[j]] == "gaussian") {
        pred.m1.e[[j]] = rnorm(nref, mean = mu.m1.e[[j]], sd = sqrt(disp.m1[[j]]))
        pred.m1.c[[j]] = rnorm(nref, mean = mu.m1.c[[j]], sd = sqrt(disp.m1[[j]]))
      }
      if (family.m1[[j]] == "poisson") {
        pred.m1.e[[j]] = rpois(nref, lambda = mu.m1.e[[j]])
        pred.m1.c[[j]] = rpois(nref, lambda = mu.m1.c[[j]])
      }
      if (family.m1[[j]] == "negbin") {
        pred.m1.e[[j]] = rnegbin(nref, mu = mu.m1.e[[j]], theta = disp.m1[[j]])
        pred.m1.c[[j]] = rnegbin(nref, mu = mu.m1.c[[j]], theta = disp.m1[[j]])
      }
      if (family.m1[[j]] == "binomial") {
        pred.m1.e[[j]] = rbinom(nref, size = 1, prob = mu.m1.e[[j]]) # size = 1 is for Bernoulli
        pred.m1.c[[j]] = rbinom(nref, size = 1, prob = mu.m1.c[[j]]) # size = 1 is for Bernoulli
      }
    }
    
    xcrep[[1]] = 1   # first exp value will be ident - use 1
    DISC[1,1] = 0
    
    ##############################
    #####   Loop of paths   ######
    ##############################
    
    c1=1
    for (c1 in 1:(npath+1)) {
      c = c1 - 1       # indexes paths
      
      # = num of rhos for noniden exp potential outcomes
      # c1 = 1 is first expected potential outcome (all exposed)
      
      m1chk = "Null"   # Init - will check for m1 in covs for each m2
      m2chk = "Null"
      
      if (c1 == 1) {m1d = 0}
      # do after 1st exp poten outcome
      if (c1 > 1) {
                
        indexDMV = matrix(indexDM[,c1-1],2,1)
        
        PrevDM = DM       # keep last DM 
        DM[indexDMV[1],indexDMV[2]] = 0   # Set one value equal to 0
        NextDM = DM        
        
        if (c1 < npath+1) {
           nindexDMV = matrix(indexDM[,c1],2,1)
           NextDM[nindexDMV[1],nindexDMV[2]] = 0
        }
        
        m1d = indexDMV[1] - 1
        m2d = indexDMV[2] - 1
        
        if (m1d > 0)  {
          m1chk = m1[[m1d]]
        } else {
          m1chk = expos[1]
        }
        if (m2d > 0) {
          m2chk = m2.eachcov[[m2d]]
        } else {
          m2chk = y.allcov
        }
        
        # Check if M1 discrep - if so use full xcorrv vector
        for (jj in 2:(numberofm1+1)) {
          for (kk in 2:(numberofm2+1)) {
            if (m1[[jj-1]] %in% m2.eachcov[[kk-1]] & DM[jj,kk] != DM[jj,1]) {  # discrep case
              DISC[1,c1] = 1          
            }
          }
        }
        
        if (DISC[1,c1] == 1) {
          xcrep[[c1]] = nxcorr
        }        
        else {xcrep[[c1]]=1}
      }  # end of C1 > 0 
      
      jdrep = xcrep[[c1]]
      
      ###################################################################
      ########  Check criterion (based on current assumption and EY   ###
      ###         to see if need to do computations                   ###
      ###################################################################
      
      #### redoing above to select paths to be plotted
      ###  just need to do for 1 bootstrap sample (use first - b=1)  
      
      ## note second condition: BOTH ds not equal - exactly or approx cancel out ##
      if (b==1) { 
        
       if (c1 > 1) {
          pathpr[[c]] = 1   #default
          
        if (i1==1 & i2<numberofm1+1){
          if ((DM[1,1]==DM[i2+1,1] & PrevDM[1,1]==PrevDM[i2+1,1]) |
              (DM[1,1]!=DM[i2+1,1] & PrevDM[1,1]!=PrevDM[i2+1,1]) )
          {pathpr[[c]] = 0}
        }
        if (i1==1 & i2>numberofm1){
          if ((DM[1,1]==DM[1,(i2-numberofm1+1)] & PrevDM[1,1]==PrevDM[1,(i2-numberofm1+1)]) |
              (DM[1,1]!=DM[1,(i2-numberofm1+1)] & PrevDM[1,1]!=PrevDM[1,(i2-numberofm1+1)]))
          {pathpr[[c]] = 0}
        }
        if (i1>1){
          if ((DM[1,i1]==DM[i2+1,i1] & PrevDM[1,i1]==PrevDM[i2+1,i1]) |
              (DM[1,i1]!=DM[i2+1,i1] & PrevDM[1,i1]!=PrevDM[i2+1,i1]))
          {pathpr[[c]] = 0}
        }
       }    
      }
          
      ##### Recanting witness sensitivity analysis loop ##################
      
      for (jd in 1:jdrep) {             
        
        if (c1 == 1 | m1chk %in% m2chk) {   #if cond not met - a null path (keep same EY)
          
          if (m1d > 0) {                   # a path through a 1st stage mediator
            nonnullm1[[m1d]][1,m2d+1] = 1    # 1 if not null  
            xcorr = xcorrv[jd]}              
          
          # Pred values for 1st stage mediators
          
          for (j in 1:numberofm1){

            # mean, pred value for M1j(dj0) (direct effects)
            j1 = j + 1
            mu.m1[[j]] = DM[j1,1]*mu.m1.e[[j]] + (1-DM[j1,1])*mu.m1.c[[j]]
            pred.m1[[j]] = DM[j1,1]*pred.m1.e[[j]] + (1-DM[j1,1])*pred.m1.c[[j]]
            
            # add next few lines for general correction, g ne 1
            
            if (i1 > 1) {
              pred.m1c[[j]] = DM[i2+1,i1]*pred.m1.e[[j]] + (1-DM[i2+1,i1])*pred.m1.c[[j]]
            }
            
            # if discrepant use
            
            mu.m12[[j]] = (1-DM[j1,1])*mu.m1.e[[j]] + DM[j1,1]*mu.m1.c[[j]]
            pred.m12[[j]] = (1-DM[j1,1])*pred.m1.e[[j]] + DM[j1,1]*pred.m1.c[[j]]
            
            # for correction term - for assumps involving an m1:
            
            mu.m1.d00[[j]] = DM[1,1] * mu.m1.e[[j]] + (1-DM[1,1]) * mu.m1.c[[j]] 

            pred.dat.y[,m1[[j]]] = pred.m1[[j]] 
          }
          
          # Pred values for 2nd stage mediators
          
          for (k in 1:numberofm2) {
            k1 = k + 1
      
            # Check whether expos is a pred for M2k
            
            if (expos[1] %in% m2.eachcov[[k]]) {
              pred.dat.m2[[k]][,expos[1]] = DM[1,k1]*exposed + (1-DM[1,k1])*control
          
              # use to get m2s for exp = d00 
              pred.dat.m2o[[k]][,expos[1]] = DM[1,1]*exposed + (1-DM[1,1])*control
              
            }
            
            # Get predicted values of 1st stage mediators for each
            #   2nd stage mediator (M2k), depending on the ds in DM
            
            for (j in 1:numberofm1) {
              
              # Check that 1st stage med is a pred for M2k
              
              if (m1[[j]] %in% m2.eachcov[[k]]) {
                
                j1 = j + 1  
                
                m1j = pred.m1[[j]]    # pred val for m1j (dir)

                if (DM[j1,k1] == DM[j1,1]) {   # non discrep case
                  pred.dat.m2[[k]][, m1[[j]]] = m1j
                  # NEW: 
                  pred.dat.m2o[[k]][, m1[[j]]] = m1j
                  
                  #### discrepant case - nonident path  ####
                  
                } else  {
                  if (first[j,jd]==1) {
                    if (family.m1[[j]] == "gaussian") {
                      ZM1 = (m1j - mu.m1[[j]])/sqrt(disp.m1[[j]])
                      ErrM1 = rnorm(nref, mean = 0, sd = 1)
                      ZM12 = (xcorr * ZM1) + sqrt(1-xcorr**2)*ErrM1
                      pred.m12[[j]] = (sqrt(disp.m1[[j]])*ZM12) + mu.m12[[j]] 
                      
                    } else {            # M1j not normal dist
                      if (family.m1[[j]] == "binomial") {
                        probjl = pbinom(m1j-1,1,mu.m1[[j]])
                        probj = pbinom(m1j, 1, mu.m1[[j]])
                      } else if (family.m1[[j]]=="poisson") {
                        probjl = ppois(m1j-1, mu.m1[[j]])
                        probj = ppois(m1j, mu.m1[[j]])
                      } else if (family.m1[[j]]== "negbin") {
                        probjl = pnbinom(m1j-1, size=disp.m1[[j]], mu=mu.m1[[j]])
                        probj = pnbinom(m1j, size=disp.m1[[j]], mu=mu.m1[[j]])
                      }
                      
                      # generate new m1j (for discr M1 cases)
                      UM1 = runif(nref, min=0, max=1)
                      probM1 = probjl + UM1 * (probj-probjl)
                      ZM1 = VGAM::probitlink(probM1)  #inv std norm
                      ErrM1 = rnorm(nref, mean = 0, sd = 1)
                      ZM12 = (xcorr * ZM1) + sqrt(1-xcorr**2)*ErrM1
                      probM12 = pnorm(ZM12)  #p(Z<ZM12)
                      
                      jj = 1 
                      prob.val = c(0)
                      
                      for (ind in 1:length(probM12)){
                        if (family.m1[[j]] == "binomial"){
                          PredM12V[[j]][ind] = as.numeric(probM12[ind] > pbinom(0,1,mu.m12[[j]][ind]))
                        } else if (family.m1[[j]] == "poisson") {
                           prob.val = 0      
                           while (prob.val < 1-tolx) {   
                            prob.val = c(prob.val, ppois(jj-1, mu.m12[[j]][ind]))
                            jj = jj + 1
                           }
                         prob.val = unique(prob.val)
                          PredM12V[[j]][ind] = as.numeric(findInterval(probM12[ind], prob.val)) - 1
                        } else if (family.m1[[j]] == "negbin") {
                          prob.val = 0      
                           while (prob.val < 1-tolx) {  
                            jj = jj + 1
                            prob.val = unique(c(prob.val, pnbinom(jj-1, size = disp.m1[[j]], mu =mu.m12[[j]][ind])))
                          }
                       PredM12V[[j]][ind] = as.numeric(findInterval(probM12[ind], unique(prob.val))) - 1
                        }
                      }
                      pred.m12[[j]] = PredM12V[[j]]
                    }
                    
                  }  # first for j, xcorr
                  
                  pred.dat.m2[[k]][, m1[[j]]] = pred.m12[[j]]
                   
                  pred.dat.m2o[[k]][, m1[[j]]] = pred.m12[[j]]
                  
                } #end of discrep case
                
              }   # if m1j a predictor of m2k
              
            }   # end of j loop (1st stage mediators)
            
            # pred values for M2k (given pred m1s, covariate values
            
            # Correction of mean of m2
            
            ####### For M2 as outcome in assumption ######
            
            mu.m2[[k]] = predict(modelb.m2[[k]], type="response",newdata=pred.dat.m2[[k]])  
           
            # m2 means for expos = d00 (need for Y vs M2 no confound assum)
            
             mu.m2o[[k]] = predict(modelb.m2[[k]], type="response",newdata=pred.dat.m2o[[k]])
             
            if (family.m2[[k]] == "gaussian") {
              pred.m2[[k]] = rnorm(nref,mean=mu.m2[[k]], sd=sqrt(disp.m2[[k]]))
            }
             
            if (family.m2[[k]] == "poisson") {
              pred.m2[[k]] = rpois(nref,lambda=mu.m2[[k]])
            }
            if (family.m2[[k]] == "binomial") {
              pred.m2[[k]] = rbinom(nref,1,prob=mu.m2[[k]])
            }
            if (family.m2[[k]] == "negbin") {
              pred.m2[[k]] = rnbinom(nref,size=disp.m2[[k]], mu=mu.m2[[k]])
            }
            
          }       # end of k loop (2nd stage mediators)
          
          #####  Expected values for Y(D)  #####
          
          pred.dat.y[,expos[1]] = DM[1,1]*exposed + (1-DM[1,1])*control
          
          for (j in 1:numberofm1) {
            pred.dat.y[,m1[[j]]] = pred.m1[[j]]
          }
          for (k in 1:numberofm2) {
            pred.dat.y[,m2[[k]]] = pred.m2[[k]]
          }
          
          # Get expected value for Y(D) and path effects (ED)
          
          # Add correction if ds (for given relationship) different
          #    for given potential outcome
       
          # Get difference in expected vals for M1 (e.g., M1(d00) vs M1(d10)) depend on path
          
          i2a = i2 - numberofm1
          
          if (i1 > 1) {
             pred.m2u = runif(nref)
          }
          
          ################################################################
          ########## Sensitivity param loop here FOR ALL ASSUMPTIONS #####
          ################################################################
          
          rhosi_vec <- numeric(nrcnt)
          
          for (rcnt in 1:nrcnt) {
                
              ri = rcnt
              rhosi = -1 + ri * deltsi
              rhosi_vec[rcnt] <- rhosi
               
           ###### For assumption involving Y  ########
            
           if (i1 == 1) {                  
              if (i2 < numberofm1 + 1) {    # assum involving 1st stage med
                m1c = pred.m1[[i2]]       
              
                if (family.m1[[i2]] != 'binomial') {
                  crct = (rhosi/stdv.m1[[i2]]) * (1/sqrt(1-rhosi**2)) *
                    (m1c*(g-1) + mu.m1.d00[[i2]] - (g*mu.m1[[i2]])) 
                }                

            if (family.y != 'binomial' & family.m1[[i2]] == 'binomial') {
              crct =  (rhosi/stdv.m1[[i2]]) * (0.5 + (1 + (dnorm(qnorm(pm1)) / stdv.m1[[i2]]))/4) * 
                 (m1c*(g-1) + mu.m1.d00[[i2]] - (g*mu.m1[[i2]]))
             }
             
             if (family.y == 'binomial' & family.m1[[i2]] == 'binomial') {
               crct =  (rhosi/stdv.m1[[i2]]) * ((dnorm(qnorm(pm1)) / stdv.m1[[i2]])**2) * 
                 (m1c*(g-1) + mu.m1.d00[[i2]] - (g*mu.m1[[i2]]))
             }
                
          } # end of assum involving 1st stage med      
             
              if (i2 > numberofm1) {       # assump involving 2nd stage med
                m2c = pred.m2[[i2a]]       
 
                if (family.m2[[i2a]] != 'binomial') {
                  crct = (rhosi/stdv.m2[[i2a]]) * (1/sqrt(1-rhosi**2)) * 
                    (m2c*(g-1) + mu.m2o[[i2a]] - (g*mu.m2[[i2a]]))        
                }
                
                if (family.y != 'binomial' & family.m2[[i2a]] == 'binomial') {
                  crct = (rhosi/stdv.m2[[i2a]]) * (0.5 + (1 + (dnorm(qnorm(pm2)) / stdv.m2[[i2a]]))/4) * 
                    (m2c*(g-1) + mu.m2o[[i2a]] - (g*mu.m2[[i2a]]))  
                }
                if (family.y == 'binomial' & family.m2[[i2a]] == 'binomial') {
                  crct = (rhosi/stdv.m2[[i2a]]) * (dnorm(qnorm(pm2)) / stdv.m2[[i2a]]) * 
                    (m2c*(g-1) + mu.m2o[[i2a]] - (g*mu.m2[[i2a]]))  
                }
              }  # end of assum involving 2nd stage med
           }    # end of assum involving Y
                
          ###### For assumption involving M2 as outcome  #######
          # if assumption does not involve Y,
          # no direct correction of mean Y - need correction for M2
                
           if (i1 > 1) {
                crct = 0
                i1m2 = i1 - 1    
                i21 = i2 + 1
                
                mu.m1.d0k[[i2]] = DM[1,i1] * mu.m1.e[[i2]] + (1-DM[1,i1]) * mu.m1.c[[i2]]
                mu.m1.djk[[i2]] = DM[i21,i1] * mu.m1.e[[i2]] + (1-DM[i21,i1]) * mu.m1.c[[i2]]
                  
                m2c = pred.m1c[[i2]]
       
                if (family.m1[[i2]] != 'binomial') {
                  crctm2 = (rhosi/stdv.m1[[i2]]) * (1/sqrt(1-rhosi**2)) *
                    (m2c*(g-1) + mu.m1.d0k[[i2]] - (g*mu.m1.djk[[i2]]))
                }
                
                if (family.y != 'binomial' & family.m1[[i2]] == 'binomial') {
                  crctm2 =  (rhosi/stdv.m1[[i2]]) * (0.5 + (1 + (dnorm(qnorm(pm1)) / stdv.m1[[i2]]))/4) * 
                    (m2c*(g-1) + mu.m1.d0k[[i2]] - (g*mu.m1.djk[[i2]]))
                }
                
                if (family.y == 'binomial' & family.m1[[i2]] == 'binomial') {
                  crctm2 =  (rhosi/stdv.m1[[i2]]) * (dnorm(qnorm(pm1)) / stdv.m1[[i2]]) * 
                    (m2c*(g-1) + mu.m1.d0k[[i2]] - (g*mu.m1.djk[[i2]]))
                }
                
                if (family.m2[[i1m2]] == "gaussian") {
                  mu.m2c[[i1m2]] = predict(modelb.m2[[i1m2]], type="response",newdata=pred.dat.m2[[i1m2]])
                  mu.m2[[i1m2]] = mu.m2c[[i1m2]] + stdv.m2[[i1m2]]*crctm2 
                  
                  pred.m2[[i1m2]] = qnorm(pred.m2u,mean=mu.m2[[i1m2]], sd=stdv.m2[[i1m2]])  
                }
          
                ###  Expected values for M2 - used if mean M2 not outcome in assump ###
                
                if (family.m2[[i1m2]] == "poisson" | family.m2[[i1m2]] == "negbin") {    
                  
                  # initiate
                  j1 = 1
                  PMj = rep(0,nref)
                  chkcum = rep(0,nref)
                  
                  # use nref - number of subjects - because evaluating simultaneously
                  # for all subjects
                  
                  while (any(pred.m2u > PMj))  {
                    j = j1 - 1  
                    j1 = j1 + 1
                    if (family.m2[[i1m2]] == "poisson") {
                      PM0j = ppois(j,mu.m2[[i1m2]])
                    }
                    else if (family.m2[[i1m2]] == "negbin") {
                      PM0j = pnbinom(j, mu=mu.m2[[i1m2]], size = disp.m2[[i1m2]]) 
                    }

                    PMj = pnorm (qnorm(PM0j) - crctm2 / sqrt(1-rhosi**2)) 
                    chkn = (pred.m2u > PMj)    # 1 when u exceeded by current PMj
                    chkcum = chkcum + chkn    # will add 1 each time u is > next cum prob
                  }
                  pred.m2[[i1m2]]=chkcum
                }
                
                if (family.m2[[i1m2]] == "binomial") {
                  PM0j = mu.m2[[i1m2]]
                  PMj = 1 - pnorm (qnorm(1 - PM0j) - crctm2 / sqrt(1-rhosi**2))
                  pred.m2[[i1m2]] = 1 * (pred.m2u > PMj)   # get binary response
                }
                
                pred.dat.y[,m2[[i1m2]]] = pred.m2[[i1m2]]         
                
           }    # if i1 > 1
                
          pred.dat.y[,m2[[k]]] = pred.m2[[k]]     
          
          ######  For ALL Assumptions  ######
                
          if (family.y == "gaussian") {    
             muY = predict(modelb.y, type="response", newdata = pred.dat.y) + stdv.y * crct
          }
                
          if (family.y == "poisson" | family.y == "negbin") {    
            muY0 = predict(modelb.y, type="response", newdata = pred.dat.y)   
            # initiate
            j1 = 0
            Pj = matrix(0,1,nref)
            sumPj = matrix(0,1,nref)
            PjT = matrix(1-tol,1,nref)   # tolerance should be close to 1
            chk = 0
            
            # use nref - number of subjects - because evaluating simultaneously
            # for all subjects
                  
            while (chk < nref & j < maxcnt)  {
              j = j1 - 1  
              j1 = j1 + 1
              if (family.y == "poisson") {
                P0j = ppois(j,muY0)
              }
              if (family.y == "negbin") {
                P0j = pnbinom(j,mu=muY0, size=disp.y)
              }
              sumPj = sumPj + Pj
              Pj = pnorm (qnorm(P0j) - crct / sqrt(1-rhosi**2)) 
              chk = sum(Pj > PjT)
            }
              
              maxY = j
              muY = (maxY * Pj) - sumPj
          }  
          
          if (family.y == "binomial") {
            muY0 = predict(modelb.y, type="response", newdata = pred.dat.y)   
            muY = 1 - pnorm (qnorm(1-muY0) - crct)
          }
          
            EYDM[[c1]][b,ri,jd] = mean(muY)   
            
          }  ######  end of sens param loop #######
          
        }    # check if non-null path
        
          else {            
               for (rcnt in 1:nrcnt) {
                 ri = rcnt  
                EYDM[[c1]][b,ri,jd] = EYDM[[c1-1]][b,ri,jd]  #null path
                }  
        nullv[1,c]= 1   # indicates path is null
        xcrep[[c1]] = 1  # just need one line   
        }
        
      }       # end of multiple xcorrs - x-world corr - (if discrep) loop
      
        # extend arrays for paths not affected by x-world correl
        
        if ((xcrep[[c1]] < nxcorr)) {     # that is - xcrep = 1, then extend
              for (lt in 2:nxcorr) {
                if (c1==1) {EYDM[[c1]][b, ,lt] = EYDM[[c1]][b, ,1]}    # EY(1)
                 else {
                  if (nullv[1,c] == 0) {EYDM[[c1]][b, ,lt] = EYDM[[c1]][b, ,1]}  # nonnull case
                   else {EYDM[[c1]][b, ,lt] = EYDM[[c1-1]][b, ,lt]}    # null case
                }
              }
         }
      
      # Get path effects as differences (or ratios) in exp values
        
     for (jd in 1:nxcorr) {      
      if (c1 > 1) {

          EDM[[c]][b,,jd] = EYDM[[c1-1]][b,,jd] - EYDM[[c1]][b,,jd]

        if (scale == "MeanDiff" || scale == "") {
          EDM_SC[[c]][b,,jd] = EYDM[[c1-1]][b,,jd] - EYDM[[c1]][b,,jd]  #pos dif dim
        }
         
        if (scale == "OR") {
          if (EYDM[[c1-1]][b,,jd] > 0 & EYDM[[c1]][b,,jd] > 0 & EYDM[[c1-1]][b,,jd] < 1 &
              EYDM[[c1]][b,,jd] < 1) {
          EDM_SC[[c]][b,,jd] = log((EYDM[[c1-1]][b,,jd]/(1-EYDM[[c1-1]][b,,jd]))
                                    /(EYDM[[c1]][b,,jd]/(1-EYDM[[c1]][b,,jd]))) # Odds Ratio option
          }
          else {stop("Expected value for Y not in (0,1) - consider a different scale")}
        }
        
        if (scale == "RR") {
          if (EYDM[[c1-1]][b,,jd] > 0 & EYDM[[c1]][b,,jd] > 0)
          EDM_SC[[c]][b,,jd] = log(EYDM[[c1-1]][b,,jd]/EYDM[[c1]][b,,jd]) # Risk Ratio option
          else {stop("Expected value for Y is negative  - use scale = 'MeanDiff'")}
        }
        
        IDENT[[c]] = 0
        if (m1d == 0)  {
          IDENT[[c]] = 1
        } else if (m1d>0) {
          if (nullv[1,c]==1 | sum(nonnullm1[[m1d]]) <= 1 ) {
            IDENT[[c]]=1}
        }
        # use IDENT to determine if need number of estimates = nxcorr
        if (IDENT[[c]] == 1) {
          xcrepd[[c]] = 1
        } else {
          xcrepd[[c]]=nxcorr
        }
        
        EDMM1[c,,b] = EDM_SC[[c]][b,,1]      
        EDMM2[c,,b,jd] = EDM_SC[[c]][b,,jd]            
        
      }     # if c1 > 1
     }      # for jd
      
    # for total

    if (scale == "MeanDiff" || scale == "") {
      EDT_SC[b,] = (EYDM[[1]][b,,1] - EYDM[[npath+1]][b,,1])  #pos dif dim total
    }

    if (scale == "OR") {
      EDT_SC[b,] = log((EYDM[[1]][b,,1]/(1-EYDM[[1]][b,,1]))/(EYDM[[npath+1]][b,,1]/
                        (1-EYDM[[npath+1]][b,,1])))    #odds ratio total
    }

    if (scale == "RR") {
      EDT_SC[b,] = log(EYDM[[1]][b,,1]/EYDM[[npath+1]][b,,1])    #risk ratio total
    }
    
    ##Define different totals using these terms
    
    # Get sums for each m1 mediator
    
    }    # end of exposure configs - paths 
    
    # Get sums for each m1 mediator (rows = bootstraps, columns = sens par)

    for (j in 1:numberofm1) {
      distinguish1 = sapply(1:ncol(EDMM1),function(x) EDMM1[,x,b] * (indexDM[1,]==j+1))
      EDM1_SC[[j]][b,] = sapply(1:ncol(distinguish1),function(x) sum(distinguish1[,x],na.rm=TRUE))  # sum for each m1, boots
    }
    
    # Get sums for each m2 mediator
        for (k in 1:numberofm2) {
          for (jd in 1:nxcorr) {
            distinguish2 = sapply(1:ncol(EDMM2),function(x) EDMM2[,x,b,jd] * (indexDM[2,]==k+1))
            EDM2_SC[[k]][b,,jd] = sapply(1:ncol(distinguish2),function(x) sum(distinguish2[,x],na.rm=TRUE))
          }
        }
    
  }    # end of bootstraps
  
#######  End of code providing estimates #######
  
######################################################################
########Start of the Creation of Plots for Sensitivity Analysis####### 
######################################################################

# Get plots for all paths

# Identify No Unmeasured Confounding Assumption
  
  if (i1 == 1) {
    av1 = y
    } else {
    av1=m2[[i1-1]]
  }
  if (i2 <= numberofm1) {
    av2 = m1[[i2]]
    } else {
    av2 = m2[[i2-numberofm1]]
  }

  # Create an empty list to store plots and tables 
plot_list <- list()
table_list <- list()

# Loop through paths and correlations
plot_index <- 1

for (pathnum in 1:npath) {  
  if (nullv[[pathnum]]==0 & pathpr[[pathnum]]==1) {   
    # output nonnull paths possibly affected by rho
  
    # get names of variables in path
    
    pathv = indexDM[,pathnum]
    m1vname=m1aug[[pathv[1]]]
    m2vname=m2aug[[pathv[2]]]
    
    # determine if xcorr relevant
    
    nxcp = nxcorr
    if (xcrep[[pathnum+1]]==1 & xcrep[pathnum]==1) {
      nxcp = 1  
    }
    
    # loop for cross corr
    
    for (xci in 1:nxcp) {
     
    EDMP <- EDM_SC[[pathnum]][ , , xci]   

  EDMP_DF <- data.frame(EDMP)

  EDMP_1_head <- tail(EDMP_DF, 1)  #last row of matrix - est based on data
  EDMP_est_SC <- data.frame(t(EDMP_1_head))   

  rho_dft <- seq(from = -1+deltsi, to = 1-deltsi, by = deltsi)
  rho_dft <- data.frame(rho_dft)

##Set up Confidence Intervals##
  CIALPHA = 1 - (CIPCT/100)    # User inputs CI percent CIPCT, e.g., 95
  CILPR = CIALPHA/2
  CIUPR = 1 - CILPR

  EDMP_boot <- t(EDMP_DF) %>%
    data.frame()
  
  EDMP_boot_q_SC <- data.frame(apply(EDMP_boot, 1, quantile, c(CILPR, CIUPR))) 
  EDMP_boot_q_SC <- t(EDMP_boot_q_SC)
  
  if (scale == "MeanDiff" || scale == "") {
    EDMP_est = EDMP_est_SC
  }
  if (scale == "OR" || scale == "RR") {
    EDMP_est = exp(EDMP_est_SC)
  }
  
  if (scale == "MeanDiff" || scale == "") {
    EDMP_boot_quant = EDMP_boot_q_SC
  }
  if (scale == "OR" || scale == "RR") {
    EDMP_boot_quant = exp(EDMP_boot_q_SC)
  }

  total_df <- cbind(rho_dft, EDMP_est, EDMP_boot_quant)

##Create the plots

  estimate_col <- paste0("X", boots + 1) 

  if (nxcp==1) {
  
  # Helper function to trim whitespace
trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x) # Remove leading and trailing whitespace
}

# Trim whitespace from variable names
m1vname <- trim_ws(m1vname)
m2vname <- trim_ws(m2vname)

# Extract the variable name from expos, ignoring the value (1)
expos_var <- expos[1]  # Only take the variable name, ignoring the numeric value

# Construct title dynamically, ensuring no extra arrow separator
title_parts <- c(expos_var)
if (nzchar(m1vname)) title_parts <- c(title_parts, m1vname)
if (nzchar(m2vname)) title_parts <- c(title_parts, m2vname)

# If path includes model.m1 (with or without m2), show xcorr
if (nzchar(m1vname)) {
  title_parts <- c(title_parts, sprintf("%s, xcorr = %s", y, xcorrv[[xci]]))
} else {
  # Only model.m2 in path: show y without xcorr
  title_parts <- c(title_parts, y)
}

title_text <- paste("Path:", paste(title_parts, collapse = " \u2192 "))
    
 plot <-  ggplot(total_df, aes_string(x = "rho_dft", y = estimate_col)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_ribbon(aes(ymin = get(paste0(round(CILPR * 100, 2), "%")),ymax = get(paste0(round(CIUPR * 100, 2), "%"))), alpha = 0.2) +
    scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                     labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                     limits = c(-1.0, 1.0)) +
      labs(title = title_text,
           x = expression(rho),     # not sure what this should be JA 6-12-25
           y = "Estimated PSE") +
        theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Adjust font size & center title
    axis.title = element_text(size = 10, face = "bold"),  # Adjust axis labels size
    axis.text = element_text(size = 10, face = "bold")    # Adjust axis tick text size
  )
 
      }
    
    else {
      
      # Helper function to trim whitespace
trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x) # Remove leading and trailing whitespace
}

# Trim whitespace from variable names
m1vname <- trim_ws(m1vname)
m2vname <- trim_ws(m2vname)

# Extract the variable name from expos, ignoring the value (1)
expos_var <- expos[1]  # Only take the variable name, ignoring the numeric value

# Construct title dynamically, ensuring no extra arrow separator
title_parts <- c(expos_var)

if (nzchar(m1vname)) title_parts <- c(title_parts, m1vname)
if (nzchar(m2vname)) title_parts <- c(title_parts, m2vname)

# If path includes model.m1 (with or without m2), show xcorr
if (nzchar(m1vname)) {
  title_parts <- c(title_parts, sprintf("%s, xcorr = %s", y, xcorrv[[xci]]))
} else {
  # Only model.m2 in path: show y without xcorr
  title_parts <- c(title_parts, y)
}

title_text <- paste("Path:", paste(title_parts, collapse = " \u2192 "))
      
      plot <-  ggplot(total_df, aes_string(x = "rho_dft", y = estimate_col)) +
        geom_line(color = "blue", linewidth = 1) +
        geom_ribbon(aes(ymin = get(paste0(round(CILPR * 100, 2), "%")),ymax = get(paste0(round(CIUPR * 100, 2), "%"))), alpha = 0.2) +
        scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                           labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                           limits = c(-1.0, 1.0)) +
      labs(title = title_text,
               x = expression(rho),
               y = 'Estimated PSE') +
          theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Adjust font size & center title
    axis.title = element_text(size = 10, face = "bold"),  # Adjust axis labels size
    axis.text = element_text(size = 10, face = "bold")    # Adjust axis tick text size
  )

       }  
 
# Store plot in the list
      plot_list[[plot_index]] <- plot
      plot_index <- plot_index + 1
 
   
  # Construct CI column names
estimate_col <- paste0("X", boots + 1)
lower_CI_col <- paste0(round(CILPR * 100, 2), "%")
upper_CI_col <- paste0(round(CIUPR * 100, 2), "%")

estimate_values   <- round(total_df[[estimate_col]], ndig)
lower_CI_values   <- round(total_df[[lower_CI_col]], ndig)
upper_CI_values   <- round(total_df[[upper_CI_col]], ndig)

# Combine into result table
result_df <- data.frame(
  Path      = title_text,
  "\u03c1"     = round(rhosi_vec, 2),
  Estimate  = estimate_values,
  `Lower CI` = lower_CI_values,
  `Upper CI` = upper_CI_values
)
 
 table_list[[length(table_list) + 1]] <- result_df  # Store table in list
 

   }  # xcorr parameter loop
  }  # determine if path to be plotted
  
 }  # end of generate plots for specific path effects
  
draw_plot_path <- function(show_titles = FALSE, start_page = 1) {
  
  # Remove NULL plots
  plot_list_clean <- Filter(Negate(is.null), plot_list)
  
  plots_per_page <- n_rows * n_cols
  n_total_plots <- length(plot_list_clean)
  
  # Start from supplied page number
  page_number <- start_page
  
  if (n_total_plots > 0) {
    
    for (i in seq(1, n_total_plots, by = plots_per_page)) {
      
      plot_subset <- plot_list_clean[
        i:min(i + plots_per_page - 1, n_total_plots)
      ]
      
      # Console title
      title_text <- paste0(
        "Estimated Path-Specific Effect (PSE) vs Sensitivity Parameter (rho)\n",
        "Assumption: No Unmeasured Confounding between '",
        av1, "' and '", av2, "'"
      )
      
      if (show_titles) {
        
        cat("\n----------------------------------------\n")
        cat("PDF Page ", page_number, "\n", sep = "")
        cat(title_text, "\n")
        
        for (p in plot_subset) {
          
          if (inherits(p, "ggplot") &&
              !is.null(p$labels$title)) {
            
            cat("  ", p$labels$title, "\n", sep = "")
          }
        }
        
        cat("----------------------------------------\n")
      }
      
      # PDF title
      title_line1 <- grid::textGrob(
        expression(
          "Estimated Path-Specific Effect (PSE) vs Sensitivity Parameter (" * rho * ")"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      title_line2 <- grid::textGrob(
        paste0(
          "Assumption: No Unmeasured Confounding between '",
          av1, "' and '", av2, "'"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      title_grob <- gridExtra::arrangeGrob(
        title_line1,
        title_line2,
        ncol = 1
      )
      
      subtitle_text <- paste0(
        "Solid line: estimates; Shaded gray area: ",
        CIPCT,
        "% confidence interval"
      )
      
      subtitle_grob <- grid::textGrob(
        subtitle_text,
        x = 0.5,
        just = "center",
        gp = grid::gpar(fontsize = 8)
      )
      
      gridExtra::grid.arrange(
        title_grob,
        subtitle_grob,
        do.call(
          gridExtra::arrangeGrob,
          c(
            plot_subset,
            ncol = n_cols,
            nrow = n_rows
          )
        ),
        ncol = 1,
        heights = c(0.12, 0.05, 0.83)
      )
      
      page_number <- page_number + 1
    }
  }
  
  # Return next available page number
  return(page_number)
}

### **Prepare Path-Specific Output Table**

# Filter out NULL or empty tables
valid_tables <- Filter(
  function(tbl) {
    !is.null(tbl) &&
      is.data.frame(tbl) &&
      nrow(tbl) > 0 &&
      any(complete.cases(tbl))
  },
  table_list
)

# Combine tables and clean Path column
if (length(valid_tables) > 0) {

  df <- dplyr::bind_rows(valid_tables) %>%
    dplyr::mutate(
      Path = stringr::str_remove(
        Path,
        "^Path:\\s*"
      ),

      # Extract xcorr from the end of Path
      Xcorr = as.numeric(
        stringr::str_extract(
          Path,
          "[-+]?[0-9]*\\.?[0-9]+$"
        )
      ),

      # Remove xcorr text from Path
      Path = stringr::str_replace(
        Path,
        ",? xcorr = [-+]?[0-9]*\\.?[0-9]+$",
        ""
      )
    ) %>%
    dplyr::select(
      Path,
      Xcorr,
      dplyr::everything()
    )

} else {

  df <- NULL
}

  ### Generate plots for 1st stage med total effects
 
# Create an empty list to store plots and tables
plot_list1 <- list()
table_list1 <- list() 

# Loop through paths and correlations
plot_index1 <- 1  
  
for (j in 1:numberofm1) {
  
  # Skip this mediator if it is not one of the assumvars
  if (!("ALL" %in% toupper(assumvars)) && !(m1[[j]] %in% assumvars)) {
    next
  } ## Added by Fang 03-19-26
    
    EDM1T <- EDM1_SC[[j]][ , ]
    
    EDM1T_DF <- data.frame(EDM1T)
    
    EDM1T_1_head <- tail(EDM1T_DF, 1)  #last row of matrix - est based on data
    EDM1T_est <- data.frame(t(EDM1T_1_head))
    
    EDM1T_boot <- t(EDM1T_DF) %>%
      data.frame()
    
    EDM1T_boot_quant <- data.frame(apply(EDM1T_boot, 1, quantile, c(CILPR, CIUPR))) 
    EDM1T_boot_quant <- t(EDM1T_boot_quant)
    total_df1 <- cbind(rho_dft, EDM1T_est, EDM1T_boot_quant)
    
    plot1 <-  ggplot(total_df1, aes_string(x = "rho_dft", y = estimate_col)) +
      geom_line(color = "blue", linewidth = 1) +
      geom_ribbon(aes(ymin = get(paste0(round(CILPR * 100, 2), "%")),ymax = get(paste0(round(CIUPR * 100, 2), "%"))), alpha = 0.2) +
      scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                         labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                         limits = c(-1.0, 1.0)) +
      labs(title = paste0("Mediator = ", m1[[j]], ""),
           x = expression(rho),
           y = 'Estimated Mediator Effect') +
        theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Adjust font size & center title
    axis.title = element_text(size = 10, face = "bold"),  # Adjust axis labels size
    axis.text = element_text(size = 10, face = "bold")    # Adjust axis tick text size
  )

    # Store plot in the list
      plot_list1[[plot_index1]] <- plot1
      plot_index1 <- plot_index1 + 1
       
# Construct CI column names
estimate_col <- paste0("X", boots + 1)
lower_CI_col <- paste0(round(CILPR * 100, 2), "%")
upper_CI_col <- paste0(round(CIUPR * 100, 2), "%")

estimate_values   <- round(total_df1[[estimate_col]], ndig)
lower_CI_values   <- round(total_df1[[lower_CI_col]], ndig)
upper_CI_values   <- round(total_df1[[upper_CI_col]], ndig)

result_df1 <- data.frame(
  Mediator   = m1[[j]],
  "\u03c1"     = round(rhosi_vec, 2),
  Estimate   = estimate_values,
  `Lower CI` = lower_CI_values,
  `Upper CI` = upper_CI_values
)


table_list1[[length(table_list1) + 1]] <- result_df1  # Store table in list

    
  }     # end of 1st stage med totals

### **Plots for 1st stage mediators**

draw_plot1 <- function(show_titles = FALSE, start_page = 1) {
  
  plots_per_page <- n_rows * n_cols
  n_total_plots <- length(plot_list1)
  
  # Start from the page number passed into the function
  page_number <- start_page
  
  if (n_total_plots > 0) {
    
    for (i in seq(1, n_total_plots, by = plots_per_page)) {
      
      # Extract plots for current page
      plot_subset1 <- plot_list1[
        i:min(i + plots_per_page - 1, n_total_plots)
      ]
      
      # Text version for Console
      title_text1 <- paste0(
        "Estimated Total Mediation Effect vs Sensitivity Parameter (rho)\n",
        "Assumption: No Unmeasured Confounding between '",
        av1, "' and '", av2, "'"
      )
      
      # Show titles in Console
      if (show_titles) {
        
        cat("\n----------------------------------------\n")
        cat("PDF Page ", page_number, "\n", sep = "")
        cat(title_text1, "\n")
        
        # Show individual plot titles
        for (p in plot_subset1) {
          
          if (inherits(p, "ggplot") &&
              !is.null(p$labels$title)) {
            
            cat("  ", p$labels$title, "\n", sep = "")
          }
        }
        
        cat("----------------------------------------\n")
      }
      
      # Main PDF title with mathematical rho
      title_line1 <- grid::textGrob(
        expression(
          "Estimated Total Mediation Effect vs Sensitivity Parameter (" * rho * ")"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      # Second title line
      title_line2 <- grid::textGrob(
        paste0(
          "Assumption: No Unmeasured Confounding between '",
          av1, "' and '", av2, "'"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      # Combine title lines
      title_grob1 <- gridExtra::arrangeGrob(
        title_line1,
        title_line2,
        ncol = 1
      )
      
      # Subtitle
      subtitle_text1 <- paste0(
        "Solid line: estimates; Shaded gray area: ",
        CIPCT,
        "% confidence interval"
      )
      
      subtitle_grob1 <- grid::textGrob(
        subtitle_text1,
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 8,
          fontface = "plain"
        )
      )
      
      # Draw current PDF page
      gridExtra::grid.arrange(
        title_grob1,
        subtitle_grob1,
        do.call(
          gridExtra::arrangeGrob,
          c(
            plot_subset1,
            ncol = n_cols,
            nrow = n_rows
          )
        ),
        ncol = 1,
        heights = c(0.12, 0.05, 0.83)
      )
      
      # Move to next page
      page_number <- page_number + 1
    }
  }
  
  # Return next available page number
  return(page_number)
}

### **Prepare First-Stage Mediator Output Table**

# Filter out NULL or empty tables
valid_tables1 <- Filter(
  function(tbl) {
    !is.null(tbl) &&
      is.data.frame(tbl) &&
      nrow(tbl) > 0
  },
  table_list1
)

if (length(valid_tables1) > 0) {
  
  combined_table1 <- dplyr::bind_rows(valid_tables1) |>
    dplyr::mutate(
      Mediator = as.character(Mediator)
    )
  
  first_stage_mediator <- unique(combined_table1$Mediator)
  
} else {
  
  combined_table1 <- NULL
  first_stage_mediator <- NULL
}

  ### Get plots for 2nd stage med totals
 # Create an empty list to store plots and tables
plot_list2 <- list()
table_list2 <- list()

# Loop through paths and correlations
plot_index2 <- 1  
  
for (k in 1:numberofm2) {
  
  # Skip this mediator if it is not one of the assumvars
  if (!("ALL" %in% toupper(assumvars)) && !(m2[[k]] %in% assumvars)) {
    next
  } 
    
    for (xci in 1:nxcorr)  {
      
      EDM2T <- EDM2_SC[[k]][ , , xci]
      EDM2T_DF <- data.frame(EDM2T)
      EDM2T_1_head <- tail(EDM2T_DF, 1)  #last row of matrix - est based on data
      EDM2T_est <- data.frame(t(EDM2T_1_head))
      
      ##Set up Confidence Intervals##

        EDM2T_boot <- t(EDM2T_DF) %>%
        data.frame() 
      
      EDM2T_boot_quant <- data.frame(apply(EDM2T_boot, 1, quantile, c(CILPR, CIUPR)))
      EDM2T_boot_quant <- t(EDM2T_boot_quant)
      total_df2 <- cbind(rho_dft, EDM2T_est, EDM2T_boot_quant)
      
      ##Create the plots
      
      title_text2 <- paste0("Mediator = ", m2[[k]], ", xcorr = ", xcorrv[[xci]], "")
      
      plot2 <-  ggplot(total_df2, aes_string(x = "rho_dft", y = estimate_col)) +
        geom_line(color = "blue", linewidth = 1) +
        geom_ribbon(aes(ymin = get(paste0(round(CILPR * 100, 2), "%")),ymax = get(paste0(round(CIUPR * 100, 2), "%"))), alpha = 0.2) +
        scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                           labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                           limits = c(-1.0, 1.0)) +
        labs(title = title_text2,
                x = expression(rho), y = 'Estimated Mediator Effect') +
               theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Adjust font size & center title
    axis.title = element_text(size = 10, face = "bold"),  # Adjust axis labels size
    axis.text = element_text(size = 10, face = "bold")    # Adjust axis tick text size
  )

      # Store plot in the list
      plot_list2[[plot_index2]] <- plot2
      plot_index2 <- plot_index2 + 1
      
      # Construct CI column names
estimate_col <- paste0("X", boots + 1)
lower_CI_col <- paste0(round(CILPR * 100, 2), "%")
upper_CI_col <- paste0(round(CIUPR * 100, 2), "%")

estimate_values   <- round(total_df2[[estimate_col]], ndig)
lower_CI_values   <- round(total_df2[[lower_CI_col]], ndig)
upper_CI_values   <- round(total_df2[[upper_CI_col]], ndig)


result_df2 <- data.frame(
  Mediator   = rep(m2[[k]], length(estimate_values)),
   Xcorr = rep(
    as.numeric(as.character(xcorrv[[xci]])),
    length(estimate_values)
  ),
  "\u03c1"     = round(rhosi_vec, 2),
  Estimate   = estimate_values,
  `Lower CI` = lower_CI_values,
  `Upper CI` = upper_CI_values
)


table_list2[[length(table_list2) + 1]] <- result_df2  # Store table in list

    }
  }     # end of 2nd stage med totals

##-- Function to draw Plot 2 pages
draw_plot2 <- function(show_titles = FALSE, start_page = 1) {
  
  plots_per_page <- n_rows * n_cols
  n_total_plots <- length(plot_list2)
  
  # Start from the page number passed into the function
  page_number <- start_page
  
  if (n_total_plots > 0) {
    
    for (i in seq(1, n_total_plots, by = plots_per_page)) {
      
      # Extract plots for current page
      plot_subset2 <- plot_list2[
        i:min(i + plots_per_page - 1, n_total_plots)
      ]
      
      # Text version for Console
      title_text2 <- paste0(
        "Estimated Total Mediation Effect vs Sensitivity Parameter (rho)\n",
        "Assumption: No Unmeasured Confounding between '",
        av1, "' and '", av2, "'"
      )
      
      # Show titles in Console
      if (show_titles) {
        
        cat("\n----------------------------------------\n")
        cat("PDF Page ", page_number, "\n", sep = "")
        cat(title_text2, "\n")
        
        # Show individual plot titles
        for (p in plot_subset2) {
          
          if (inherits(p, "ggplot") &&
              !is.null(p$labels$title)) {
            
            cat("  ", p$labels$title, "\n", sep = "")
          }
        }
        
        cat("----------------------------------------\n")
      }
      
      # Main PDF title with mathematical rho
      title_line1 <- grid::textGrob(
        expression(
          "Estimated Total Mediation Effect vs Sensitivity Parameter (" * rho * ")"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      # Second title line
      title_line2 <- grid::textGrob(
        paste0(
          "Assumption: No Unmeasured Confounding between '",
          av1, "' and '", av2, "'"
        ),
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        )
      )
      
      # Combine title lines
      title_grob2 <- gridExtra::arrangeGrob(
        title_line1,
        title_line2,
        ncol = 1
      )
      
      # Subtitle
      subtitle_text2 <- paste0(
        "Solid line: estimates; Shaded gray area: ",
        CIPCT,
        "% confidence interval"
      )
      
      subtitle_grob2 <- grid::textGrob(
        subtitle_text2,
        x = 0.5,
        just = "center",
        gp = grid::gpar(
          fontsize = 8,
          fontface = "plain"
        )
      )
      
      # Draw current PDF page
      gridExtra::grid.arrange(
        title_grob2,
        subtitle_grob2,
        do.call(
          gridExtra::arrangeGrob,
          c(
            plot_subset2,
            ncol = n_cols,
            nrow = n_rows
          )
        ),
        ncol = 1,
        heights = c(0.12, 0.05, 0.83)
      )
      
      # Move to next page
      page_number <- page_number + 1
    }
  }
  
  # Return next available page number
  return(page_number)
}

# Combine Tables
valid_tables2 <- Filter(
  function(tbl) {
    !is.null(tbl) &&
      nrow(tbl) > 0 &&
      any(complete.cases(tbl))
  },
  table_list2
)

# Combine tables
if (length(valid_tables2) > 0) {
  
  df2 <- dplyr::bind_rows(valid_tables2)
  
  # Clean mediator labels
  df2 <- df2 %>%
    dplyr::mutate(
      Mediator = stringr::str_remove(
        as.character(Mediator),
        "^Mediator\\s*=\\s*"
      ),
      Xcorr = as.numeric(as.character(Xcorr))
    ) %>%
    dplyr::select(
      Mediator,
      Xcorr,
      dplyr::everything()
    )
  
} else {
  
  df2 <- NULL
}

##-- Collect this pair's tables and draw this pair's plots, so that every pair
##-- ends up in the same CSV files and in the same PDF
pair_label <- paste0(av1, " & ", av2)

if (is.data.frame(df) && nrow(df) > 0) {
  df$Rho.Between <- pair_label
  all_df[[length(all_df) + 1]] <- df
}
if (is.data.frame(combined_table1) && nrow(combined_table1) > 0) {
  combined_table1$Rho.Between <- pair_label
  all_df1[[length(all_df1) + 1]] <- combined_table1
}
if (is.data.frame(df2) && nrow(df2) > 0) {
  df2$Rho.Between <- pair_label
  all_df2[[length(all_df2) + 1]] <- df2
}

if (!is.null(pdf_dev)) {
  grDevices::dev.set(pdf_dev)
  page_number <- draw_plot_path(show_titles = TRUE, start_page = page_number)
  page_number <- draw_plot1(show_titles = TRUE, start_page = page_number)
  page_number <- draw_plot2(show_titles = TRUE, start_page = page_number)
}

# end of second variable in confounded pair
    }  
  }   # end of first variable in confounded pair
  
# Sensitivity Analysis Output
cat("----------------------------------------\n")
cat("Exposure ('exposed' value):",paste0(expos[1],"(", expos[2], ")", "\n"))
cat("Final Outcome:", y, "\n")
cat("Sample Size Used:", nrow(gdata), "\n")
cat("Number of Bootstrap Samples / Number Used:", sprintf("%d / %d\n", nrow(EDT) - 1, boots))

if (!identical(ref, "entire sample")) {
  cat("Reference Group:", ref[1], "= (", paste(ref[-1], collapse = ", "), "), ",
      "Reference Group Multiplier (refmult):", refmult, "\n")
} else {
  cat("Reference Group: Whole Sample, refmult:", refmult, "\n")
}
cat("Scale:", scale, "\n")
if (!is.null(cluster)) {
  cat("Cluster:", cluster, "\n")
}

cat("Confidence Interval Percentage:", CIPCT, "\n")
cat("g (between-correlation factor):", g, "\n")
cat("----------------------------------------\n")

##-- Export the output tables (all assumption pairs combined)
df_out  <- if (length(all_df)  > 0) dplyr::bind_rows(all_df)  else NULL
df1_out <- if (length(all_df1) > 0) dplyr::bind_rows(all_df1) else NULL
df2_out <- if (length(all_df2) > 0) dplyr::bind_rows(all_df2) else NULL

if (table_out[1] == "true") {
  if (nchar(table_out[2]) == 0) {
    print("WARNING: Please Specify Path to Save the Output Files")
  } else {
    
    # Create the output folder if it does not exist yet
    if (!dir.exists(table_out[2])) {
      dir.create(table_out[2], recursive = TRUE, showWarnings = FALSE)
    }
    
    ### Path specific
    out_path <- glue::glue(table_out[2], "/gmedsens_pathspec_out.csv")
    if (!is.null(df_out) && nrow(df_out) > 0) {
      write.csv(df_out, out_path, row.names = FALSE)
    }
    
    ### First stage mediator
    out_path1 <- glue::glue(table_out[2], "/gmedsens_firststage_out.csv")
    if (!is.null(df1_out) && nrow(df1_out) > 0) {
      write.csv(df1_out, out_path1, row.names = FALSE)
    }
    
    ### Second stage mediator
    out_path2 <- glue::glue(table_out[2], "/gmedsens_secondstage_out.csv")
    if (!is.null(df2_out) && nrow(df2_out) > 0) {
      write.csv(df2_out, out_path2, row.names = FALSE)
    }
  }
}

##-- Close the plots PDF and confirm that it was written
if (!is.null(pdf_dev)) {
  if (pdf_dev %in% grDevices::dev.list()) grDevices::dev.off(pdf_dev)
  if (file.exists(out_plot) && file.info(out_plot)$size > 0) {
    cat("Plots saved to: ", normalizePath(out_plot), "\n", sep = "")
  } else {
    cat("WARNING: The plots PDF was not created:", out_plot, "\n")
    warning("The plots PDF was not created: ", out_plot, call. = FALSE)
  }
}

invisible(NULL)
  }  #end of gmedsens function
