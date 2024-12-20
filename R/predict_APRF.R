

#' Obtain predictions from a Asset Pricing Random Forest (APRF).
#'
#' @param object : a \code{APRF} output of train_APRF() function.
#' @param X [matrix]: data.frame of the fixed effects for the new observations to be predicted.
#' @param id [vector]: vector of the identifiers of the new observations to be predicted.
#' @param time [vector]: vector of the time measurements of the new observations to be predicted.
#' @param ... : low levels arguments.
#'
#' @import stats
#' @import ranger
#' @importFrom rpart
#' @return vector of the predicted output for the new observations.
#'
#' @export
#'
#' @examples\donttest{}
#'
predict.APRF <- function(object, X, id, time, ...) {
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  Z <- matrix(1, nrow = nrow(X), ncol = 1)
  if (class(object$forest) == "rpart") {
    f <- predict(object$forest, data = data.frame(X))
  }
  else if (class(object$forest) == "ranger") {
    f <- ranger::predictions(predict(object$forest,data = data.frame(X = X)))
  }
  Time <- object$time
  id_btilde <- object$id_btilde
  Ypred <- rep(0, length(id))
  id.app <- object$id
  if (is.null(object$sto)) return(f)
  if (object$sto == "none") {
    for (i in 1:length(unique(id))) {
      w <- which(id == unique(id)[i])
      k <- which(id_btilde == unique(id)[i])
      Ypred[w] <- f[w] + Z[w, , drop = FALSE] %*% object$random_effects[k, ]
    }
    return(Ypred)
  }

  for (i in 1:length(unique(id))) {
    w <- which(id == unique(id)[i])
    k <- which(id_btilde == unique(id)[i])
    om <- which(id.app == unique(id)[i])
    Ypred[w] <- f[w] + Z[w, , drop = FALSE] %*% object$random_effects[k, ] + predict.sto(object$omega[om], Time[om], time[w], object$sto)
  }
  return(Ypred)
}
