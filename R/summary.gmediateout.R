#' @title Gmediate Object Summaries
#'
#' @description This function generates result summary table of gmediate function output object
#' @param object output object of gmediate function
#' @param ... arguments passed to other methods
#' @export
#' @return NULL

summary.gmediateout = function(object, ...){
  ##-- Estimated path effects --##
  cat("\n\n")
  cat("Causal Mediation Analysis", "\n")
  cat("Exposure:", object$expos[1], "\n")
  cat("Outcome:", object$y, "\n")
  cat("Sample Size Used:", dim(object$data)[1], "\n")
  cat("Number of Bootstrap Samples Used:", dim(object$EDT)[1]-1, "/", object$bootsims, "\n")
  if (!is.null(object$ref)){
    cat("Reference Group:", object$ref[1], "= (", object$ref[2:length(object$ref)], "),", "Reference Group Multiplier (refmult) =", object$refmult, "\n")
  }
  if (is.null(object$ref)){
    cat("Reference Group: Whole Sample,", "Reference Group Multiplier (refmult)  =", object$refmult, "\n")
  }
  if (!is.null(object$cluster)){
    cat("Cluster:", object$cluster)
  }
  
  cat("\n\n")
  cat("-------------------------------\n")
  cat("Mediation/Path Effect Estimates\n")
  cat("-------------------------------")
  cat("\n\n")
  cat("Individual Path Effects")
  cat("\n")
  
  
  print(object$indiv.path, row.names = F)
  if (sum(object$nullv) > 0){
    cat(" * = a priori null path")
  }
  
  cat("\n\n")
  cat("Total Effect")
  cat("\n")
  cat("")
  
  print(object$total.path, row.names = F)
  cat("\n\n")
  
  cat("First Stage Mediators")
  cat("\n")
  cat("")
  
  print(object$first.path, row.names = F)
  cat("\n\n")
  
  cat("Second Stage Mediators")
  cat("\n")
  cat("")
  
  print(object$second.path, row.names = F)
  cat("\n\n")
  
  ##-- Proportion of Estimated Path Effects --##
  ##-- Draw Table 5 : Proportion of estimated path effects
  cat("\n")
  cat("----------------------------------------------\n")
  cat("Estimated Proportions of Total Effect Mediated\n")
  cat("----------------------------------------------")
  cat("\n\n")
  cat("Individual Path Effects")
  cat("\n")
  
  print(object$prop.indiv.path, row.names = F)
  if (sum(object$nullv) > 0){
    cat(" * = a priori null path")
  }
  cat("\n\n")
  
  ##-- Draw Table 6 : Proportion of estimated 1st stage mediator effects
  cat("First Stage Mediators")
  cat("\n")
  cat("")
  print(object$prop.first.path, row.names = F)
  cat("\n\n")
  
  ##-- Draw Table 7 : Proportion of estimated 2st stage mediator effects
  cat("Second Stage Mediators")
  cat("\n")
  cat("")
  
  print(object$prop.second.path, row.names = F)
  cat("\n\n")
}
