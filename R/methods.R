#' @export
predict.lmm <- function(object, new_data, ...){
  if (!inherits(object, "lmm"))
    stop("object not of class \"lmm\"")

  args <- list(...)

  if(length(args)==0){ # if the grouping variable is not known for newdata, I simply use the average estimates
    return(c(object$beta %*% t(new_data)))
  }

  grouping_variable <- unlist(new_data[, args$grouping_variable_name])
  model_matrix_new_data <-
    subset(new_data, select = -get(args$grouping_variable_name))
  random_effect_component <- object$mu_raneff[grouping_variable]

  c(object$beta %*% t(model_matrix_new_data)) +
    random_effect_component
}

#' @export
predict.mlmm <- function(object, new_data, ...){
  if (!inherits(object, "mlmm"))
    stop("object not of class \"mlmm\"")

  args <- list(...)

  if(length(args)==0){ # if the grouping variable is not known for newdata, I simply use the average estimates
    return(new_data %*%object$BETA)
  }
  grouping_variable <- unlist(new_data[, args$grouping_variable_name])

  model_matrix_new_data <-
    subset(new_data, select = -get(args$grouping_variable_name))

  random_effect_component <- t(object$mu_raneff[,,grouping_variable])

  model_matrix_new_data %*%object$BETA +
    random_effect_component
}
