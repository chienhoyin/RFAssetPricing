#' Train Longitudinal Random Forest Models for Asset Pricing
#'
#' @import ranger
#' @import rpart
#'
#'
#' @param model Character string specifying the model type. Must be one of:
#'   \itemize{
#'     \item "mert": Mixed Effects Random Tree
#'     \item "smert": Stochastic Mixed Effects Random Tree
#'     \item "merf": Mixed Effects Random Forest
#'     \item "smerf": Stochastic Mixed Effects Random Forest
#'     \item "rf": Standard Random Forest
#'     \item "reemt": Random Effects EM Tree
#'     \item "sreemt": Stochastic Random Effects EM Tree
#'   }
#' @param train_X A data.frame containing the fixed effects predictors.
#' @param train_Y Numeric vector containing the response variable.
#' @param train_time Numeric vector containing the time measurements associated with each observation.
#' @param train_id Vector containing the identifiers for different trajectories/entities.
#' @param n_iter Integer specifying the maximum number of iterations for model fitting. Default is 100.
#' @param n_trees Integer specifying the number of trees to grow (for forest-based models). Default is 500.
#' @param mtry Integer specifying the number of variables randomly sampled at each split. Default is ceiling(ncol(train_X) / 3).
#' @param n_cores Integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param sto Character string specifying the stochastic process type (for stochastic variants). Must be one of:
#'   \itemize{
#'     \item "none": No stochastic process (default)
#'     \item "BM": Brownian motion
#'   }
#' @param del Numeric value specifying the convergence threshold. The algorithm stops when the relative
#'   difference in log-likelihood between iterations is smaller than this value. Default is 0.001.
#'
#' @return For model="rf", returns a ranger object (see \link[ranger]{ranger} for details).
#' For all other models, returns a fitted model object of class "APRF" containing:
#'   \itemize{
#'     \item forest: The fitted random forest (ranger object) or tree (rpart object) depending on model type
#'     \item random_effects: Predictions of random effects for different trajectories
#'     \item var_random_effects: Estimated variance-covariance matrix of random effects
#'     \item sigma: Estimated residual variance parameter
#'     \item id_btilde: Identifiers associated with random effects predictions
#'     \item omega: (For stochastic variants) Predicted stochastic processes
#'     \item sigma_sto: (For stochastic variants) Estimated volatility parameter
#'     \item time: Time measurements vector
#'     \item sto: Stochastic process specification
#'     \item Vraisemblance: Log-likelihood values across iterations
#'     \item id: Vector of trajectory identifiers
#'     \item OOB: Out-of-bag error estimates (for forest-based models)
#'   }
#' See \link[ranger]{ranger} and \link[rpart]{rpart} for details about the underlying forest/tree objects.
#'
#' @examples
#' \dontrun{
#' # Train a MERF model
#' model <- train_APRF(
#'   model = "merf",
#'   train_X = as.matrix(train_data[, char_names]),
#'   train_Y = train_data$xret,
#'   train_time = as.numeric(as.factor(train_data$date)),
#'   train_id = train_data$gvkey,
#'   n_iter = 100,
#'   n_trees = 500
#' )
#'
#' # Train a stochastic variant with Brownian motion
#' model_sto <- train_APRF(
#'   model = "smerf",
#'   train_X = as.matrix(train_data[, char_names]),
#'   train_Y = train_data$xret,
#'   train_time = as.numeric(as.factor(train_data$date)),
#'   train_id = train_data$gvkey,
#'   sto = "BM"
#' )
#' }
#'
#' @export
train_APRF <- function(
    model = c("mert", "smert", "merf", "smerf", "rf", "reemt", "sreemt"),
    train_X,
    train_Y,
    train_time,
    train_id,
    n_iter = 100,
    n_trees = 500,
    mtry = NULL,
    n_cores = 1,
    sto = "none",
    del = 0.001) {
  # Input validation
  model <- match.arg(model)

  if (is.null(mtry)) {
    mtry <- ceiling(ncol(train_X) / 3)
  }

  # Convert X to appropriate format based on model
  if (model %in% c("mert", "reemt", "sreemt")) {
    train_X <- as.data.frame(train_X)
  } else {
    train_X <- as.matrix(train_X)
  }

  # Set number of threads for ranger
  #old_threads <- options("ranger.num.threads")[[1]]
  #on.exit(options(ranger:ranger.num.threads = old_threads))
  #options(ranger::ranger.num.threads = n_cores)

  train_Z <- matrix(1, nrow = length(train_Y), ncol = 1)

  # Train the specified model
  trained_model <- switch(model,
                          "mert" = MERT(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            sto = "none",
                            delta = del
                          ),
                          "smert" = MERT(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            sto = sto,
                            delta = del
                          ),
                          "merf" = MERF(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            mtry = mtry,
                            ntree = n_trees,
                            sto = "none",
                            n_threads = n_cores,
                            delta = del
                          ),
                          "smerf" = MERF(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            mtry = mtry,
                            ntree = n_trees,
                            sto = sto,
                            n_threads = n_cores,
                            delta = del
                          ),
                          "rf" = list(forest=ranger::ranger(
                            data = data.frame(X = train_X, y = train_Y),
                            dependent.variable.name = "y",
                            mtry = mtry,
                            num.trees = n_trees,
                            importance = "impurity",
                            oob.error = TRUE,
                            num.threads = n_cores)
                          ),
                          "reemt" = REEMtree(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            sto = "none",
                            delta = del
                          ),
                          "sreemt" = REEMtree(
                            X = train_X,
                            Y = train_Y,
                            Z = train_Z,
                            id = train_id,
                            time = train_time,
                            iter = n_iter,
                            sto = sto,
                            delta = del
                          )
  )
  class(trained_model) <- 'APRF'
  return(trained_model)
}
