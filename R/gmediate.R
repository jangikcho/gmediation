#' @title Generalized Causal Mediation and Path Analysis
#'
#' @description This function estimates causal mediation and path effects for models with two stages of mediation
#' @param data a data frame
#' @param model.m1 a list of the model fit from each first stage mediators
#' @param model.m2 a list of the model fit from each second stage mediators
#' @param model.y an object of the model fit from the response variable
#' @param expos a vector with the name of the exposure variable and exposure value e.g. c("exposure variable name", c(exposure value))
#' @param ref a vector with the name of the reference group variable and reference group value e.g. c("reference group variable name", c(reference group value))
#' @param refmult numeric value used as a multiplier to increase the number of monte carlo draws
#' @param seed seed to generate bootstrap samples
#' @param bootsims number of bootstrap samples to generate
#' @param cluster the cluster variable that will be used for cluster bootstrap resampling
#' @param decomp option with two choices for decomposition of total exposure effect
#' @param sig.level confidence level for bootstrap confidence interval
#' @param sens.par a vector with numeric values for sensitivity parameters to be used when the path is non-identifiable
#' @importFrom stats ecdf getCall na.omit pbinom pnbinom pnorm ppois predict quantile rbinom rnbinom rnorm rpois runif
#' @importFrom plyr count
#' @import MASS
#' @export
#' @seealso \href{https://cran.r-project.org/package=mediation}{mediation}
#' @return NULL
#' @details
#' This function carries out inference for mediation (path) effects based on user-specified models. The function allows for two stages of mediation (up to three links from exposure/treatment, X, to final outcome, Y), and an arbitrary number of mediators at each stage. The user specifies an allowed generalized regression model for each outcome (Y and each mediator). Allowed distributions (link functions) include: normal (identity link), Bernoulli (logit link), Poisson (log link), and negative binomial (log link); note that canonical links are used for each distribution. Unsaturated models (i.e., leaving out earlier stage mediators) are allowed, possibly imposing a prior null paths. Baseline (pre-exposure) covariates may be included in each model. The user-specified models are re-fit in the function using the complete cases (based on all the model variables).
#'
#' The user may also specify a reference group, effectively defining the empirical multivariate baseline covariate distribution used in the mediation formula; this may be viewed as defining the sub-population to which inference is performed. The default reference group is the whole (complete case) sample. The reference group may also be ‘cloned’ by a specified integer multiple (“refmult”) to reduce Monte Carlo error (at the possible cost of longer computing time).
#'
#' An extended mediation formula is used to estimate each path effect. The sum of effects (or ‘overall effect’) for paths through each mediator is also computed. For paths (and overall) effects that are ‘non-identifiable’, multiple estimates are provided according to specified sensitivity parameters values (“sens.par”, a list of values between – 1 and 1), the default values being -0.9, 0, 0.5, and 0.9. Path effects have alternative versions depending on the manner in which the total exposure effect is decomposed. Two choices, using “decomp”, are provided; the ordering of mediators in specified models also generally affect the decomposition.
#'
#' Bootstrap resampling is used to obtain bootstrap percentile method confidence intervals (at the 95 percent, or user-specified, confidence level). The option “cluster” can be used to specify a unit (aside from the individual, that is, each data record, which is the default) for which cluster bootstrap resampling is performed.
#'
#' @examples
#' library(gmediation)
#'
#' ## mediation analysis for two stage mediators with one mediator at each stage
#' model.y1 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + seal, family = binomial,
#'                 data = dental)
#' model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental)
#' model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush, family = binomial,
#'                  data = dental)
#' single <- gmediate(data = dental, model.m1 = list(model.m11), model.m2 = list(model.m21),
#'                    model.y = model.y1, expos = c("ses", 1), ref = c("ses", c(0,1)), seed = 1234,
#'                    bootsims = 100, cluster = NULL, decomp = 1, sig.level = 95, sens.par = c(0),
#'                    refmult = 1)
#' single
#' names(single)
#'
#'\dontrun{
#' ## mediation analysis for two stage mediators with two mediator at each stage
#' model.y2 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + visit + seal + ohi,
#'                 family = binomial, data = dental)
#' model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental)
#' model.m12 <- glm(visit ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental)
#' model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush + visit, family = binomial,
#'                  data = dental)
#' model.m22 <- glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + visit, family = gaussian,
#'                  data = dental)
#' gmediate(data = dental, model.m1 = list(model.m11,model.m12), model.m2 = list(model.m21,model.m22),
#'          model.y = model.y2, expos = c("ses", 1), ref = c("ses", c(0,1)), seed = 352,
#'          bootsims = 100, cluster = NULL, decomp = 2, sig.level = 95, sens.par = c(0,0.9),
#'          refmult = 2)
#'}
#'



gmediate = function(data, model.m1, model.m2, model.y, expos, ref = NULL,
                    refmult = 1, seed = 1234, bootsims = 500, cluster = NULL,
                    decomp = 1, sig.level = 95, sens.par = c(-0.9, 0, 0.5, 0.9)){

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
  ## Define Warning messages ##
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

  # number of correlations
  nrho = length(sens.par)

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
  pred.dat.m1.c = list()
  pred.dat.m1.e = list()

  ##-- Define reference group
  if (is.null(ref)){
    data.cr = gdata
  } else {
    data.cr0 = c()
    for (i in 2:length(ref)){
      data.cr0 = rbind(data.cr0, gdata[gdata[,ref[1]] == ref[i],])
    }

    data.cr = c()
    if(refmult == 1){
      data.cr = data.cr0
    } else {
      for (i in 1:refmult){
        data.cr = rbind(data.cr0, data.cr)
      }
    }
  }

  ##-- Define dataset for the reference group in Y
  pred.dat.y = data.cr[,y.allcov]

  ##-- Define dataset for the reference group in the 2nd stage mediators
  for (k in 1:numberofm2) {
    pred.dat.m2[[k]] = data.cr[,m2.eachcov[[k]]]
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

  ####################################
  ##-- Generate Bootstrap Samples --##
  ####################################
  bootsamples = vector("list", bootsims + 1)

  if (is.null(cluster)){
    nclus = nrow(gdata)
    gdata$numsamp = c(1:nclus)

    for (i in 1:bootsims){
      samp = sample(gdata$numsamp, size = nclus, replace = TRUE)
      bootsamples[[i]] = gdata[samp,]
    }
    bootsamples[[bootsims+1]] = gdata
  }

  if (!is.null(cluster)){
    nclus = length(unique(gdata[,cluster]))

    for (i in 1:bootsims){
      samp = sample(data[,cluster], size = nclus, replace = TRUE)
      bootsamples[[i]] = data.frame()
      for (j in 1:length(samp)){
        bootsamples[[i]] = rbind(bootsamples[[i]],data[data[,cluster] == samp[j],])
      }
    }
    bootsamples[[bootsims+1]] = gdata
  }

  ##########################
  ##-- Initialize lists --##
  ##########################
  modelb.y = c()
  modelb.m2 = list()
  modelb.m1 = list()

  mu.m2 = list()
  mu.m1 = list()
  mu.m12 = list()
  mu.m1.e = list()
  mu.m1.c = list()

  pred.m2 = list()
  pred.m12 = list()
  pred.m1.e = list()
  pred.m1.c = list()
  pred.m1 = list()
  PredM12V = vector("list", numberofm1)

  disp.m2 = list()    # list with dispersion par values 2nd stage
  disp.m1 = list()    # list with dispersion par values 1st stage


  rhrep = list()
  rhrepd = list()
  nonnullm1 = list()  # added 9-9-16
  IDENT = list()      # added 9-9-16

  #initialize nonnullm1 457-459 added 9-9-16
  for (i in 1:numberofm1) {
    nonnullm1[[i]] = matrix(0,1,numberofm2+1)
  }

  EYDM = list()      # exp pot outcome matrices for boots, rhos
  EDM = list()       # expected path effects
  EDMM2 = vector("list", nrho)
  EDV = list()       # first element (rho) for each path effect
  EST = list()
  EDLB = vector("list", npath)      # CI lower bounds
  EDUB = vector("list", npath)      # CI upper bounds
  PV = vector("list", npath)       # P values
  EDB1 = list()
  EDB2 = vector("list", numberofm2)
  for (t in 1:numberofm2){
    EDB2[[t]] = matrix(0, nrho, bootsims+1)
  }

  EDT = matrix(0, bootsims+1, 1)

  rho = c()

  EDM1 = list()
  EDM1L = list()
  EDM1U = list()
  PV1 = list()

  EDM2 = vector("list", numberofm2)
  EDM2L = vector("list", numberofm2)
  EDM2U = vector("list", numberofm2)
  PV2 = vector("list", numberofm2)

  ########################
  ##-- Bootstrap loop --##
  ########################
  #initialize for each EY
  for (p in 1:(npath+1)) {
    EYDM[[p]] = matrix(0,bootsims+1,nrho)
  }

  for (c in 1 : npath){
    EDM[[c]] = matrix(rep(0,((bootsims+1)*nrho)),bootsims+1,nrho)  #init
  }

  for (b in 1:(bootsims + 1)){

    nullv = matrix(0,1,npath)
    DISC = matrix(0,1,npath+1)

    ##--  Fit models for each sample (boot and orig)
    for (k in 1:numberofm2) {
      call.m2[[k]]$data = bootsamples[[b]][, c(m2[[k]], m2.eachcov[[k]])]

      if (family.m2[[k]] == "negbin"){
        modelb.m2[[k]] = glm.nb(call.m2[[k]]$formula, call.m2[[k]]$data)
      } else {
        modelb.m2[[k]] = eval.parent(call.m2[[k]])
      }

      if (family.m2[[k]] == "gaussian" | family.m2[[k]] == "negbin") {
        disp.m2[[k]] = summary(modelb.m2[[k]])$disp
      }
      # Get means for indivs in ref grp
      mu.m2[[k]] = predict(modelb.m2[[k]], type = "response", newdata = pred.dat.m2[[k]])
    }

    for (j in 1:numberofm1){
      # Get data (all model vars) for bootstrap sample for refit

      # Note: indi (used below) not yet defined
      # should be the individual indices acc to cluster boot draw
      call.m1[[j]]$data = bootsamples[[b]][, c(m1[[j]], m1.eachcov[[j]])]

      # Refit model using current sample (boot or orig)
      if (family.m1[[j]] == "negbin"){
        modelb.m1[[j]] = glm.nb(call.m1[[j]]$formula, call.m1[[j]]$data)
      } else {
        modelb.m1[[j]] = eval.parent(call.m1[[j]])
      }

      if (family.m1[[j]] == "gaussian" | family.m1[[j]] == "negbin") {
        disp.m1[[j]] = summary(modelb.m1[[j]])$disp
      }

      # Get means for indivs in ref grp under each expos status
      mu.m1.e[[j]] = predict(modelb.m1[[j]], type = "response", newdata = pred.dat.m1.e[[j]])
      mu.m1.c[[j]] = predict(modelb.m1[[j]], type = "response", newdata = pred.dat.m1.c[[j]])
    }

    call.y$data = bootsamples[[b]][, c(y, y.allcov)]
    if (family.y == "negbin"){
      modelb.y = glm.nb(call.y$formula, call.y$data)
    } else {
      modelb.y = eval.parent(call.y)    # refit Y model
    }


    ##-- End of fitting models

    ###################################################################
    ###    Med Formula: Draw values for each med, each path eff     ###
    ###################################################################

    DM = matrix(1, numberofm1+1, numberofm2+1)

    first = matrix(1, numberofm1, nrho)

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

    rhrep[[1]] = 1   # first exp value will be ident - use 1
    DISC[1,1] = 0

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
        DM[indexDMV[1],indexDMV[2]] = 0   # Set one value equal to 0

        # Remove from here - put to 602-614

        # indic mediator (each stage) involved in current exp poten out
        m1d = indexDMV[1] - 1
        m2d = indexDMV[2] - 1

        # fixed the following 6 lines 8-10-16
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

        # Check if M1 discrep - if so use full rho vector
        for (jj in 2:(numberofm1+1)) {
          for (kk in 2:(numberofm2+1)) {
            if (m1[[jj-1]] %in% m2.eachcov[[kk-1]] & DM[jj,kk] != DM[jj,1]) {   # discrep case
              DISC[1,c1] = 1          # added 8-10-16
            }
          }
        }

        if (DISC[1,c1] == 1) {
          rhrep[[c1]] = nrho
        }        # added 8-11-16
        else {rhrep[[c1]]=1}
      }

      jdrep = rhrep[[c1]]

      for (jd in 1:jdrep) {

        if (c1 == 1 | m1chk %in% m2chk) {   #if cond not met - a null path (keep same EY)

          if (m1d > 0) {                   # a path through a 1st stage mediator
            nonnullm1[[m1d]][1,m2d+1] = 1    # 1 if not null  made correction 9-9-16
            rho = sens.par[jd]}

          # Pred values for 1st stage mediators

          for (j in 1:numberofm1){

            # mean, pred value for M1j(dj0) (direct effects)
            j1 = j + 1
            mu.m1[[j]] = DM[j1,1]*mu.m1.e[[j]] + (1-DM[j1,1])*mu.m1.c[[j]]
            pred.m1[[j]] = DM[j1,1]*pred.m1.e[[j]] + (1-DM[j1,1])*pred.m1.c[[j]]

            # if discrepant use - add next line 11-1-16
            mu.m12[[j]] = (1-DM[j1,1])*mu.m1.e[[j]] + DM[j1,1]*mu.m1.c[[j]]
            pred.m12[[j]] = (1-DM[j1,1])*pred.m1.e[[j]] + DM[j1,1]*pred.m1.c[[j]]
          }

          # Pred values for 2nd stage mediators

          for (k in 1:numberofm2) {
            k1 = k + 1

            # Check whether expos is a pred for M2k

            if (expos[1] %in% m2.eachcov[[k]]) {
              pred.dat.m2[[k]][,expos[1]] = DM[1,k1]*exposed + (1-DM[1,k1])*control
            }

            # Get predicted values of 1st stage mediators for each
            #   2nd stage mediator (M2k), depending on the ds in DM

            for (j in 1:numberofm1) {

              # Check that 1st stage med is a pred for M2k

              if (m1[[j]] %in% m2.eachcov[[k]]) {

                j1 = j + 1    # uncommented 11-8-16

                m1j = pred.m1[[j]]    # pred val for m1j (dir)

                if (DM[j1,k1] == DM[j1,1]) {   # non discrep case - fix 11-1-16
                  pred.dat.m2[[k]][, m1[[j]]] = m1j

                  #### discrepant case - nonident path  ####

                } else  {
                  if (first[j,jd]==1) {
                    if (family.m1[[j]] == "gaussian") {
                      ZM1 = (m1j - mu.m1[[j]])/sqrt(disp.m1[[j]])
                      ErrM1 = rnorm(nref, mean = 0, sd = 1)
                      ZM12 = (rho * ZM1) + sqrt(1-rho**2)*ErrM1
                      pred.m12[[j]] = (sqrt(disp.m1[[j]])*ZM12) + mu.m12[[j]]  #edited 11-1-16

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
                      ZM1 = VGAM::probit(probM1)  #inv std norm
                      ErrM1 = rnorm(nref, mean = 0, sd = 1)
                      ZM12 = (rho * ZM1) + sqrt(1-rho**2)*ErrM1
                      probM12 = pnorm(ZM12)  #p(Z<ZM12)

                      for (ind in 1:length(probM12)){
                        if (family.m1[[j]] == "binomial"){
                          PredM12V[[j]][ind] = as.numeric(probM12[ind] > pbinom(0,1,mu.m12[[j]][ind]))
                        } else if (family.m1[[j]] == "poisson") {
                          prob.val = c()
                          for (jj in 1:100) {
                            prob.val = c(prob.val, ppois(jj-1, mu.m12[[j]][ind]))
                          }
                          prob.val = unique(prob.val)
                          PredM12V[[j]][ind] = as.numeric(findInterval(probM12[ind], prob.val)) - 1
                        } else if (family.m1[[j]] == "negbin") {
                          prob.val = c()
                          for (jj in 1:100) {
                            prob.val = unique(c(prob.val, pnbinom(jj-1, size = disp.m1[[j]], mu =mu.m12[[j]][ind])))
                          }
                          PredM12V[[j]][ind] = as.numeric(findInterval(probM12[ind], unique(prob.val))) - 1
                        }
                      }
                      pred.m12[[j]] = PredM12V[[j]]
                    }

                  }  # first for j, rho

                  pred.dat.m2[[k]][, m1[[j]]] = pred.m12[[j]]

                } #end of discrep case

              }   # if m1j a predictor of m2k

            }   # end of j loop (1st stage mediators)

            # pred values for M2k (given pred m1s, covariate values

            mu.m2[[k]] = predict(modelb.m2[[k]], type="response",newdata=pred.dat.m2[[k]])
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
          muY = predict(modelb.y, type="response", newdata = pred.dat.y)

          EYDM[[c1]][b,jd] = mean(muY)

        }

        else {EYDM[[c1]][b,1] = EYDM[[c1-1]][b,1]  #null path  added 8-11-16
        nullv[1,c]= 1   # indicates path is null  added 8-10-16
        rhrep[[c1]] = 1  # just need one line     added 8-11-16
        }

      }       # end of multiple rhos (if discrep) loop

      # modified following paragraph 8-11-16

      if ((rhrep[[c1]] < nrho)) {     # that is - rhrep = 1, then extend
        for (lt in 2:nrho) {
          if (c1==1) {EYDM[[c1]][b,lt] = EYDM[[c1]][b,1]}    # EY(1)
          else {
            if (nullv[1,c] == 0) {EYDM[[c1]][b,lt] = EYDM[[c1]][b,1]}  # nonnull case
            else {EYDM[[c1]][b,lt] = EYDM[[c1-1]][b,lt]}    # null case
          }
        }
      }

      # Get path effects as differences in exp values
      if (c1 > 1) {
        if (b==1) {
          if (c==1) {
            EDMM1 = matrix(rep(0,(npath*(bootsims+1))),npath,bootsims+1)  #init
          }
        }
        EDM[[c]][b,] = EYDM[[c1-1]][b,] - EYDM[[c1]][b,]  #pos dif dim
        #rhrepd[[c]] = max(rhrep[[c1]],rhrep[[c1-1]])  deleted this 9-9-16

        # added following (lines 801-813) 9-9-16
        IDENT[[c]] = 0
        if (m1d == 0)  {
          IDENT[[c]] = 1
        } else if (m1d>0) {
          if (nullv[1,c]==1 | sum(nonnullm1[[m1d]]) <= 1 ) {
            IDENT[[c]]=1}
        }
        # use IDENT to determine if need number of estimates = nrho
        if (IDENT[[c]] == 1) {
          rhrepd[[c]] = 1
        } else {
          rhrepd[[c]]=nrho
        }

        EDMM1[c,b] = EDM[[c]][b,1]      # turn into matrix - first rho

        #        if (rhrepd[c]==1) {
        #          EDM[[c]][b,] = EDM[[c]][b,] %*% rep(1,nrho)  # expand
        #        }

        for (jd in 1:nrho) {
          if (b==1 & c==1) {
            EDMM2[[jd]] = matrix(rep(0,(npath*(bootsims+1))),npath,bootsims+1)  #init
          }
          EDMM2[[jd]][c,b] = EDM[[c]][b,jd]
        }
      }
    }           # end of exposure configs

    # for total
    EDT[b,] = (EYDM[[1]][b,] - EYDM[[npath+1]][b,])[rhrep[[1]]]

    # Get sums for each m1 mediator
    for (j in 1:numberofm1) {
      distinguish1 = sapply(1:ncol(EDMM1),function(x) EDMM1[,x] * (indexDM[1,]==j+1))
      EDB1[[j]] = sapply(1:ncol(distinguish1),function(x) sum(distinguish1[,x],na.rm=TRUE))  # sum for each m1, boots
    }

    # Get sums for each m2 mediator
    for (k in 1:numberofm2) {
      for (jd in 1:nrho) {
        distinguish2 = sapply(1:ncol(EDMM2[[jd]]),function(x) EDMM2[[jd]][,x] * (indexDM[2,]==k+1))
        EDB2[[k]][jd,] = sapply(1:ncol(distinguish2),function(x) sum(distinguish2[,x],na.rm=TRUE))
      }
    }
  }    # end of bootstraps

  CIALPHA = 1 - (sig.level/100)    # User inputs CI percent, e.g., 95
  CILPR = CIALPHA/2
  CIUPR = 1 - CILPR

  # Need m1, m2 name vectors to be extended with 1st entry for direct effect
  m1m = matrix(data = m1,1,numberofm1)     # m1names is list of m1 names
  m2m = matrix(data = m2,1,numberofm2)     # m2names            m2
  rhom = matrix(data = sens.par,1,nrho)
  m1ex = cbind("   ", m1m)
  m2ex = cbind("   ", m2m)
  rhoex = cbind("   ", rhom)   # Need user to specify vector rho for m1 vars

  #########################################
  ##-- Organize Estimated Path Effects --##
  #########################################
  for (c in 1:npath) {
    EST[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])  #init for each
    EDLB[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])  #init for each
    EDUB[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])  #init for each
    PV[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])  #init for each

    for (jd in 1:rhrepd[[c]]) {
      EST[[c]][1,jd] = round(EDM[[c]][bootsims+1,jd],4)   # Est path effect orig data
      EDLB[[c]][1,jd] = round(quantile(EDM[[c]][1:bootsims,jd],CILPR,na.rm = TRUE),4) # CI lb
      EDUB[[c]][1,jd] = round(quantile(EDM[[c]][1:bootsims,jd],CIUPR,na.rm = TRUE),4) # CI ub

      # get p-value for test of effect = 0
      PV[[c]][1,jd] = format(round((2 * min(ecdf(EDM[[c]][1:bootsims,jd])(0), ecdf(-EDM[[c]][1:bootsims,jd])(0))),4),nsmall = 4)
    }    # over jd (multiple rhos)
  }   # over c paths

  ESTT = round(EDT[bootsims+1,],4)   # Est total effect orig data
  ETLB = round(quantile(EDT[1:bootsims,],CILPR,na.rm = TRUE),4) # CI lb
  ETUB = round(quantile(EDT[1:bootsims,],CIUPR,na.rm = TRUE),4)  # CI ub
  PVT = format(round(2 * min(ecdf(EDT[1:bootsims,])(0), ecdf(-EDT[1:bootsims,])(0)),4),nsmall = 4)

  for (j in 1:numberofm1) {
    EDM1[[j]] = round(EDB1[[j]][bootsims+1],4)
    EDM1L[[j]] = round(quantile(EDB1[[j]][1:bootsims],CILPR,na.rm = TRUE),4)  # CI lower bound
    EDM1U[[j]] = round(quantile(EDB1[[j]][1:bootsims],CIUPR,na.rm = TRUE),4)  # CI upper bound

    # get p-value for test of effect = 0
    PV1[[j]] = format(round(2 * min(ecdf(EDB1[[j]][1:bootsims])(0), ecdf(-EDB1[[j]][1:bootsims])(0)),4),nsmall = 4)
  }

  for (k in 1:numberofm2) {
    for (jd in 1:nrho) {
      EDM2[[k]][jd] = round(EDB2[[k]][jd,bootsims+1],4)
      EDM2L[[k]][jd] = round(quantile(EDB2[[k]][jd,1:bootsims],CILPR),4)  # CI lb
      EDM2U[[k]][jd] = round(quantile(EDB2[[k]][jd,1:bootsims],CIUPR),4)  # CI ub

      # get p-value for test of effect = 0
      PV2[[k]][jd] = format(round(2 * min(ecdf(EDB2[[k]][jd,1:bootsims])(0),ecdf(-EDB2[[k]][jd,1:bootsims])(0)),4),nsmall = 4)
    }
  }

  ###############################
  ##-- Print out the results --##
  ###############################
  cat("\n\n")
  cat("Causal Mediation Analysis", "\n")
  cat("Exposure:", expos[1], "\n")
  cat("Outcome:", y, "\n")
  cat("Sample Size Used:", dim(gdata)[1], "\n")
  cat("Number of Bootstrap Samples Used:", dim(EDT)[1]-1, "/", bootsims, "\n")
  if (!is.null(ref)){
    cat("Reference Group:", ref[1], "= (", ref[2:length(ref)], "),", "Reference Group Multiplier (refmult) =", refmult, "\n")
  }
  if (is.null(ref)){
    cat("Reference Group: Whole Sample,", "Reference Group Multiplier (refmult)  =", refmult, "\n")
  }
  if (!is.null(cluster)){
    cat("Cluster:", cluster)
  }

  cat("\n\n")
  cat("-------------------------------\n")
  cat("Mediation/Path Effect Estimates\n")
  cat("-------------------------------")
  cat("\n\n")
  cat("Individual Path Effects")
  cat("\n")

  effect1p = c()
  effect2p = c()
  rhop = c()
  effectsp = c()
  effectup = c()
  effectlp = c()
  pvp = c()

  for (c in 1:npath){
    for (jd in 1:rhrepd[[c]]){
      if (nullv[c] == 1){
        EST[[c]][1,jd] = "*"
        EDLB[[c]][1,jd] = "*"
        EDUB[[c]][1,jd] = "*"
        PV[[c]][1,jd] = "*"
      }
    }
  }

  for (c in 1:npath){
    for (jd in 1:rhrepd[[c]]){
      effect1p[length(effect1p) + 1] = paste(m1ex[indexDM[1,c]])
      effect2p[length(effect2p) + 1] = paste(m2ex[indexDM[2,c]])
      if (nullv[c] == 0 & IDENT[[c]] == 1 || nullv[c] == 1){
        rhop[length(rhop) + 1] = paste("-")
      } else {
        rhop[length(rhop) + 1] = paste(rhom[jd])
      }
      effectsp[length(effectsp) + 1] = EST[[c]][1,jd]
      effectlp[length(effectlp) + 1] = EDLB[[c]][1,jd]
      effectup[length(effectup) + 1] = EDUB[[c]][1,jd]
      pvp[length(pvp) + 1] = PV[[c]][1,jd]
    }
  }

  arrow = c(rep("--->",length(effect1p)))

  table1 = data.frame(effect1p,arrow,effect2p,rhop,effectsp,effectlp,effectup,pvp)

  colnames(table1) = c(" ", "Path", " ", "Rho", "Estimate", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table1, row.names = F)
  if (sum(nullv) > 0){
    cat(" * = a priori null path")
  }

  cat("\n\n")
  cat("Total Effect")
  cat("\n")
  cat("")
  table2 = data.frame(ESTT,ETLB,ETUB,PVT)
  colnames(table2) = c("Estimate", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table2, row.names = F)
  cat("\n\n")

  cat("First Stage Mediators")
  cat("\n")
  cat("")
  m1names = c()
  m1est = c()
  m1ciu = c()
  m1cil = c()
  m1pv = c()
  for (i in 1:numberofm1){
    m1names[i] = paste(m1[[i]])
    m1est[i] = EDM1[[i]]
    m1cil[i] = EDM1L[[i]]
    m1ciu[i] = EDM1U[[i]]
    m1pv[i] = PV1[[i]]
  }
  table3 = data.frame(m1names,m1est,m1cil,m1ciu,m1pv)
  colnames(table3) = c("Mediator","Estimate", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table3, row.names = F)
  cat("\n\n")

  cat("Second Stage Mediators")
  cat("\n")
  cat("")
  m2names = c()
  m2est = c()
  m2ciu = c()
  m2cil = c()
  m2pv = c()
  for (i in 1:numberofm2){
    for (j in 1:nrho){
      m2names[length(m2names) + 1] = paste(m2[[i]])
      m2est[length(m2est) + 1] = EDM2[[i]][j]
      m2cil[length(m2cil) + 1] = EDM2L[[i]][j]
      m2ciu[length(m2ciu) + 1] = EDM2U[[i]][j]
      m2pv[length(m2pv) + 1] = PV2[[i]][j]
    }
  }
  m2rho = c(rep(sens.par,numberofm2))
  table4 = data.frame(m2names,m2rho,m2est,m2cil,m2ciu,m2pv)
  colnames(table4) = c("Mediator","Rho","Estimate", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table4, row.names = F)
  cat("\n\n")

  #######################################################
  ##-- Organize Proportion of Estimated Path Effects --##
  #######################################################
  ##-- Calculate proportion of estimated path effects from the original data
  PEST = vector("list", npath)
  for (c in 1:npath){
    PEST[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])
    for (j in 1:rhrepd[[c]]){
      PEST[[c]][j] = round(EST[[c]][j]/ESTT,4)
    }
  }

  ##-- Calculate proportions from bootstrap samples to get 95% CI
  PESTCI = vector("list", npath)
  for (c in 1:npath){
    PESTCI[[c]] = matrix(rep(0,rhrepd[[c]]),bootsims,rhrepd[[c]])
    for (k in 1:bootsims){
      for (jd in 1:rhrepd[[c]]){
        PESTCI[[c]][k,jd] = EDM[[c]][k,jd]/EDT[[k]]
      }
    }
  }

  ##-- Select 95% CI
  PESTCILB = vector("list", npath)
  PESTCIUB = vector("list", npath)
  PPV = vector("list", npath)
  for (c in 1:npath){
    PESTCILB[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])
    PESTCIUB[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])
    PPV[[c]] = matrix(rep(0,rhrepd[[c]]),1,rhrepd[[c]])
    for (jd in 1:rhrepd[[c]]){
      PESTCILB[[c]][1,jd] = round(quantile(PESTCI[[c]][1:bootsims,jd],CILPR,na.rm = TRUE),4) # CI lb
      PESTCIUB[[c]][1,jd] = round(quantile(PESTCI[[c]][1:bootsims,jd],CIUPR,na.rm = TRUE),4) # CI ub
      # get p-value for test of proportion of effect = 0
      PPV[[c]][1,jd] = format(round((2 * min(ecdf(PESTCI[[c]][1:bootsims,jd])(0), ecdf(-PESTCI[[c]][1:bootsims,jd])(0))),4),nsmall = 4)
    }
  }

  ##-- Draw Table 5 : Proportion of estimated path effects
  cat("\n")
  cat("----------------------------------------------\n")
  cat("Estimated Proportions of Total Effect Mediated\n")
  cat("----------------------------------------------")
  cat("\n\n")
  cat("Individual Path Effects")
  cat("\n")

  proeffect1p = c()
  proeffect2p = c()
  prorhop = c()
  proeffectsp = c()
  proeffectup = c()
  proeffectlp = c()
  propvp = c()

  for (c in 1:npath){
    for (jd in 1:rhrepd[[c]]){
      if (nullv[c] == 1){
        PEST[[c]][1,jd] = "*"
        PESTCILB[[c]][1,jd] = "*"
        PESTCIUB[[c]][1,jd] = "*"
        PPV[[c]][1,jd] = "*"
      }
    }
  }

  for (c in 1:npath){
    for (jd in 1:rhrepd[[c]]){
      proeffect1p[length(proeffect1p) + 1] = paste(m1ex[indexDM[1,c]])
      proeffect2p[length(proeffect2p) + 1] = paste(m2ex[indexDM[2,c]])
      if (nullv[c] == 0 & IDENT[[c]] == 1 || nullv[c] == 1){
        prorhop[length(prorhop) + 1] = paste("-")
      } else {
        prorhop[length(prorhop) + 1] = paste(rhom[jd])
      }
      proeffectsp[length(proeffectsp) + 1] = PEST[[c]][1,jd]
      proeffectlp[length(proeffectlp) + 1] = PESTCILB[[c]][1,jd]
      proeffectup[length(proeffectup) + 1] = PESTCIUB[[c]][1,jd]
      propvp[length(propvp) + 1] = PPV[[c]][1,jd]
    }
  }

  arrow = c(rep("--->",length(proeffect1p)))

  table5 = data.frame(proeffect1p,arrow,proeffect2p,prorhop,proeffectsp,proeffectlp,proeffectup,propvp)

  colnames(table5) = c(" ", "Path", " ", "Rho", "Est Prop", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table5, row.names = F)
  if (sum(nullv) > 0){
    cat(" * = a priori null path")
  }
  cat("\n\n")

  ##-- Calculate proportion of estimated effects of 1st stage mediators
  PEDM1 = vector("list",numberofm1)
  for (j in 1:numberofm1){
    PEDM1[[j]] = round(EDM1[[j]]/ESTT,4)
  }

  ##-- Calculate proportions from bootstrap samples to get 95% CI
  PEDMCI1 = vector("list",numberofm1)
  for (j in 1:numberofm1){
    for (b in 1:bootsims){
      PEDMCI1[[j]][b] = EDB1[[j]][b]/EDT[b]
    }
  }

  ##-- Select 95% CI
  PEDM1L = vector("list", numberofm1)
  PEDM1U = vector("list", numberofm1)
  PPV1 = vector("list", numberofm1)
  for (j in 1:numberofm1){
    PEDM1L[[j]] = round(quantile(PEDMCI1[[j]][1:bootsims],CILPR,na.rm = TRUE),4)  # CI lower bound
    PEDM1U[[j]] = round(quantile(PEDMCI1[[j]][1:bootsims],CIUPR,na.rm = TRUE),4)  # CI upper bound
    # get p-value for test of effect = 0
    PPV1[[j]] = format(round(2 * min(ecdf(PEDMCI1[[j]][1:bootsims])(0), ecdf(-PEDMCI1[[j]][1:bootsims])(0)),4),nsmall = 4)
  }

  ##-- Draw Table 6 : Proportion of estimated 1st stage mediator effects
  cat("First Stage Mediators")
  cat("\n")
  cat("")
  prom1names = c()
  prom1est = c()
  prom1ciu = c()
  prom1cil = c()
  prom1pv = c()

  for (i in 1:numberofm1){
    prom1names[i] = paste(m1[[i]])
    prom1est[i] = PEDM1[[i]]
    prom1cil[i] = PEDM1L[[i]]
    prom1ciu[i] = PEDM1U[[i]]
    prom1pv[i] = PPV1[[i]]
  }

  table6 = data.frame(prom1names,prom1est,prom1cil,prom1ciu,prom1pv)
  colnames(table6) = c("Mediator","Est Prop", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table6, row.names = F)
  cat("\n\n")

  ##-- Calculate proportion of estimated effects of 2nd stage mediators
  PEDM2 = vector("list",numberofm2)
  for (k in 1:numberofm2){
    for (jd in 1:nrho){
      PEDM2[[k]][jd] = round(EDM2[[k]][jd]/ESTT,4)
    }
  }

  ##-- Calculate proportions from bootstrap samples to get 95% CI
  PEDMCI2 = vector("list",numberofm2)
  for (k in 1:numberofm2){
    PEDMCI2[[k]] = matrix(rep(0,bootsims),nrho,bootsims)
    for (jd in 1:nrho){
      for (b in 1:bootsims){
        PEDMCI2[[k]][jd,b] = EDB2[[k]][jd,b]/EDT[b]
      }
    }
  }

  ##-- Select 95% CI
  PEDM2L = vector("list", numberofm2)
  PEDM2U = vector("list", numberofm2)
  PPV2 = vector("list", numberofm2)
  for (k in 1:numberofm2){
    for (jd in 1:nrho){
      PEDM2L[[k]][jd] = round(quantile(PEDMCI2[[k]][jd,1:bootsims],CILPR,na.rm = TRUE),4)  # CI lower bound
      PEDM2U[[k]][jd] = round(quantile(PEDMCI2[[k]][jd,1:bootsims],CIUPR,na.rm = TRUE),4)  # CI upper bound
      # get p-value for test of effect = 0
      PPV2[[k]][jd] = format(round(2 * min(ecdf(PEDMCI2[[k]][jd,1:bootsims])(0), ecdf(-PEDMCI2[[k]][jd,1:bootsims])(0)),4),nsmall = 4)
    }
  }

  ##-- Draw Table 7 : Proportion of estimated 2st stage mediator effects
  cat("Second Stage Mediators")
  cat("\n")
  cat("")

  prom2names = c()
  prom2est = c()
  prom2ciu = c()
  prom2cil = c()
  prom2pv = c()
  for (i in 1:numberofm2){
    for (j in 1:nrho){
      prom2names[length(prom2names) + 1] = paste(m2[[i]])
      prom2est[length(prom2est) + 1] = PEDM2[[i]][j]
      prom2cil[length(prom2cil) + 1] = PEDM2L[[i]][j]
      prom2ciu[length(prom2ciu) + 1] = PEDM2U[[i]][j]
      prom2pv[length(prom2pv) + 1] = PPV2[[i]][j]
    }
  }
  prom2rho = c(rep(sens.par,numberofm2))
  table7 = data.frame(prom2names,prom2rho,prom2est,prom2cil,prom2ciu,prom2pv)
  colnames(table7) = c("Mediator","Rho","Est Prop", paste(sig.level, "% CI Lower", sep=""), paste(sig.level, "% CI Upper", sep=""), "p-value")
  print(table7, row.names = F)
  cat("\n\n")

  ## prepare for output objects
  # each output table
  table1n = table1[,-2]
  colnames(table1n) = c("1st Stage Mediator","2nd Stage Mediators",colnames(table1)[4:8])
  table2n = table2
  table3n = table3
  table4n = table4
  table5n = table5[,-2]
  colnames(table1n) = c("1st Stage Mediator","2nd Stage Mediator",colnames(table1)[4:8])
  table6n = table6
  table7n = table7

  # bootstrap estimates
  EDM.out = vector("list",npath)
  for (c in 1:npath){
    EDM.out[[c]] = EDM[[c]][-(bootsims+1),]
  }

  # model fit
  mfit = vector("list", numberofm1+numberofm2+1)
  # model for first stage mediators
  for (j in 1:numberofm1){
    mfit[[j]][[1]] = modelb.m1[[j]]$formula
    mfit[[j]][[2]] = summary(modelb.m1[[j]])$coefficient
  }
  # model for second stage mediators
  for (k in 1:numberofm2){
    mfit[[numberofm1+k]][[1]] = modelb.m2[[k]]$formula
    mfit[[numberofm1+k]][[2]] = summary(modelb.m2[[k]])$coefficient
  }
  # model for y
  mfit[[numberofm1+numberofm2+1]][[1]] = modelb.y$formula
  mfit[[numberofm1+numberofm2+1]][[2]] = summary(modelb.y)$coefficient

  return(invisible(list(indiv.path = table1n, total.path = table2n, first.path = table3n, second.path = table4n,
              prop.indiv.path = table5n, prop.first.path = table6n, prop.second.path = table7n,
              boot.est = EDM.out, model.fit = mfit)))

}  #end of gmediate function


