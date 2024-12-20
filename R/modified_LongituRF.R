#' Title
#'
#'
#'
#' @import stats
#' @import ranger
#'
#' @keywords internal
MERF <- function(X, Y, id, Z, iter = 100, mtry = ceiling(ncol(X) / 3), ntree = 500, time, sto, delta = 0.001, n_threads = 1 ) {
  #old_threads <- options("ranger.num.threads")[[1]]
  #on.exit(options(ranger::ranger.num.threads = old_threads))
  #options(ranger::ranger.num.threads = n_threads)

  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0, nind, q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1, q)) ### init
  epsilonhat <- rep(0, length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0, length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL

  if (class(sto) == "character") {

    if (sto == "none") {
      for (i in 1:iter) {
        print(i)
        ystar <- rep(NA, length(Y))
        for (k in 1:nind) { #### on retrace les effets al?atoires
          indiv <- which(id == unique(id)[k])
          ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }

        forest <- ranger::ranger(
          data = data.frame(X = X, y = ystar),
          dependent.variable.name = "y",
          mtry = mtry,
          num.trees = ntree,
          importance = "impurity",
          oob.error = TRUE,
          num.threads = n_threads
        )
        fhat <- ranger::predictions(predict(forest, data = data.frame(X = X)))
        OOB[i] <- forest$prediction.error

        for (k in 1:nind) {
          indiv <- which(id == unique(id)[k])
          V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv))
          btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }

        sigm <- sigmahat
        sigmahat <- sig(sigma = sigmahat, id = id, Z = Z, epsilon = epsilonhat, Btilde = Btilde)
        Btilde <- bay(bhat = btilde, Bhat = Btilde, Z = Z, id = id, sigmahat = sigm)
        Vrai <- c(Vrai, logV(Y, fhat, Z, time, id, Btilde, 0, sigmahat, sto))

        if (i > 1) inc <- abs((Vrai[i - 1] - Vrai[i]) / Vrai[i - 1])
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest = forest, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, Vraisemblance = Vrai, id = id, time = time, OOB = OOB)
          class(sortie) <- "APRF"
          return(sortie)
        }
      }
      sortie <- list(forest = forest, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, Vraisemblance = Vrai, id = id, time = time, OOB = OOB)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }
  for (i in 1:iter) {
    print(i)
    ystar <- rep(0, length(Y))
    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }

    forest <- ranger::ranger(
      data = data.frame(X = X, y = ystar),
      dependent.variable.name = "y",
      mtry = mtry,
      num.trees = ntree,
      importance = "impurity",
      oob.error = TRUE,
      num.threads = 20
    )
    fhat <- ranger::predictions(predict(forest, data = data.frame(X = X)))
    OOB[i] <- forest$prediction.error

    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      K <- sto_analysis(sto, time[indiv])
      V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv)) + sigma2 * K
      btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
      omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat, id, Z, epsilonhat, Btilde, time, sigma2, sto)
    Btilde <- bay_sto(btilde, Btilde, Z, id, sigm, time, sigma2, sto)
    sigma2 <- gam_sto(sigm, id, Z, B, time, sigma2, sto, omega)
    Vrai <- c(Vrai, logV(Y, fhat, Z[, , drop = FALSE], time, id, Btilde, sigma2, sigmahat, sto))

    if (i > 1) inc <- abs((Vrai[i - 1] - Vrai[i]) / Vrai[i - 1])
    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest = forest, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id, OOB = OOB)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }
  sortie <- list(forest = forest, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id, OOB = OOB)
  class(sortie) <- "APRF"
  return(sortie)
}

#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
bay_sto <- function(bhat, Bhat, Z, id, sigmahat, time, sigma2, sto) { #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind) {
    w <- which(id == unique(id)[j])
    K <- sto_analysis(sto, time[w])
    V <- Z[w, , drop = FALSE] %*% Bhat %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigmahat), length(w), length(w)) + sigma2 * K
    D <- D + (bhat[j, ] %*% t(bhat[j, ])) + (Bhat - Bhat %*% t(Z[w, , drop = FALSE]) %*% solve(V) %*% Z[w, , drop = FALSE] %*% Bhat)
  }
  D <- D / nind
  return(D)
}

#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
sig_sto <- function(sigma, id, Z, epsilon, Btilde, time, sigma2, sto) { #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind) {
    w <- which(id == unique(id)[j])
    K <- sto_analysis(sto, time[w])
    V <- Z[w, , drop = FALSE] %*% Btilde %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigma), length(w), length(w)) + sigma2 * K
    sigm <- sigm + t(epsilon[w]) %*% epsilon[w] + sigma * (length(w) - sigma * (sum(diag(solve(V)))))
  }
  sigm <- sigm / Nombre
  return(sigm)
}

#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
gam_sto <- function(sigma, id, Z, Btilde, time, sigma2, sto, omega) {
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind) {
    indiv <- which(id == unique(id)[k])
    K <- sto_analysis(sto, time[indiv])
    V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigma), length(indiv), length(indiv)) + sigma2 * K
    Omeg <- omega[indiv]
    gam <- gam + (t(Omeg) %*% solve(K) %*% Omeg) + sigma2 * (length(indiv) - sigma2 * sum(diag(solve(V) %*% K)))
  }
  return(as.numeric(gam) / Nombre)
}


#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
predict.sto <- function(omega, time.app, time.test, sto) {
  pred <- rep(0, length(time.test))

  if (class(sto) == "function") {
    for (i in 1:length(time.test)) {
      inf <- which(time.app <= time.test[i])
      sup <- which(time.app > time.test[i])
      if (length(inf) > 0) {
        if (length(sup) > 0) {
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app == time_inf)], omega[which(time.app == time_sup)]))
        }
        if (length(sup) == 0) {
          time_inf <- max(time.app[inf])
          pred[i] <- omega[which(time.app == time_inf)] * ((sto(time.test[i], max(time.app))) / sto(max(time.app), max(time.app)))
        }
      }
      if (length(sup) > 0 & length(inf) == 0) {
        time_sup <- min(time.app[sup])
        pred[i] <- omega[which(time.app == time_sup)] * ((sto(time.test[i], min(time.app))) / sto(min(time.app), min(time.app)))
      }
    }
    return(pred)
  } else {
    for (i in 1:length(time.test)) {
      inf <- which(time.app <= time.test[i])
      sup <- which(time.app > time.test[i])
      if (length(inf) > 0) {
        if (length(sup) > 0) {
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app == time_inf)], omega[which(time.app == time_sup)]))
        }
        if (length(sup) == 0) {
          time_inf <- max(time.app[inf])
          if (sto == "BM") {
            pred[i] <- omega[which(time.app == time_inf)]
          }
          if (sto == "OrnUhl") {
            pred[i] <- omega[which(time.app == time_inf)] * (exp(-abs(time.test[i] - max(time.app)) / 2))
          }
          if (sto == "BBridge") {
            pred[i] <- omega[which(time.app == time_inf)] * ((1 - time.test[i]) / (1 - max(time.app)^2))
          }
        }
      }
      if (length(sup) > 0 & length(inf) == 0) {
        time_sup <- min(time.app[sup])
        if (sto == "BM") {
          pred[i] <- omega[which(time.app == time_sup)] * (time.test[i] / min(time.app))
        }
        if (sto == "OrnUhl") {
          pred[i] <- omega[which(time.app == time_sup)] * (exp(-abs(time.test[i] - min(time.app)) / 2))
        }
        if (sto == "BBridge") {
          pred[i] <- omega[which(time.app == time_sup)] * (time.test[i] / min(time.app))
        }
      }
    }
  }
  return(pred)
}
#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
sto_analysis <- function(sto, time) {
  MAT <- matrix(0, length(time), length(time))

  if (class(sto) == "function") {
    for (i in 1:length(time)) {
      for (j in 1:length(time)) {
        MAT[i, j] <- sto(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto == "BM") {
    for (i in 1:length(time)) {
      for (j in 1:length(time)) {
        MAT[i, j] <- min(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto == "OrnUhl") {
    for (i in 1:length(time)) {
      for (j in 1:length(time)) {
        MAT[i, j] <- exp(-abs(time[i] - time[j]) / 2)
      }
    }
    return(MAT)
  }

  if (sto == "BBridge") {
    for (i in 1:length(time)) {
      for (j in 1:length(time)) {
        MAT[i, j] <- min(time[i], time[j]) - time[i] * time[j]
      }
    }
    return(MAT)
  }
}

#' Title
#'
#'
#'
#'
#' @import stats
#'
#' @keywords internal
sig <- function(sigma, id, Z, epsilon, Btilde) { #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind) {
    w <- which(id == unique(id)[j])
    V <- Z[w, , drop = FALSE] %*% Btilde %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigma), length(w), length(w))
    sigm <- sigm + t(epsilon[w]) %*% epsilon[w] + sigma * (length(w) - sigma * (sum(diag(solve(V)))))
  }
  sigm <- sigm / Nombre
  return(sigm)
}
#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat, Bhat, Z, id, sigmahat) { #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind) {
    w <- which(id == unique(id)[j])
    V <- Z[w, , drop = FALSE] %*% Bhat %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigmahat), length(w), length(w))
    D <- D + (bhat[j, ] %*% t(bhat[j, ])) + (Bhat - Bhat %*% t(Z[w, , drop = FALSE]) %*% solve(V) %*% Z[w, , drop = FALSE] %*% Bhat)
  }
  D <- D / nind
  return(D)
}
#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
Moy <- function(id, Btilde, sigmahat, Phi, Y, Z) {
  S1 <- 0
  S2 <- 0
  nind <- length(unique(id))
  for (i in 1:nind) {
    w <- which(id == unique(id)[i])
    V <- Z[w, , drop = FALSE] %*% Btilde %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigmahat), length(w), length(w))
    S1 <- S1 + t(Phi[w, , drop = FALSE]) %*% solve(V) %*% Phi[w, , drop = FALSE]
    S2 <- S2 + t(Phi[w, , drop = FALSE]) %*% solve(V) %*% Y[w]
  }
  epsilon <- 1e-6
  dim_S1 <- dim(S1)
  S1_reg <- S1 + epsilon * diag(nrow = dim_S1, ncol = dim_S1)

  # Return regularized solution
  return(solve(S1_reg) %*% S2)
}
#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
Moy_sto <- function(id, Btilde, sigmahat, Phi, Y, Z, sto, time, sigma2) {
  S1 <- 0
  S2 <- 0
  nind <- length(unique(id))
  for (i in 1:nind) {
    w <- which(id == unique(id)[i])
    K <- sto_analysis(sto, time[w])
    V <- Z[w, , drop = FALSE] %*% Btilde %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigmahat), length(w), length(w)) + sigma2 * K
    S1 <- S1 + t(Phi[w, , drop = FALSE]) %*% solve(V) %*% Phi[w, , drop = FALSE]
    S2 <- S2 + t(Phi[w, , drop = FALSE]) %*% solve(V) %*% Y[w]
  }
  return(solve(S1) %*% S2)
}

#' Title
#'
#'
#'
#' @import stats
#' @importFrom rpart predict
#'
#' @keywords internal
MERT <- function(X, Y, id, Z, iter = 100, time, sto, delta = 0.001) {
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0, nind, q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1, q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0, length(Y))
  sigma2 <- 1
  inc <- 1
  Vrai <- NULL
  id_omega <- sto

  if (class(sto) == "character") {
    if (sto == "none") {
      for (i in 1:iter) {
        print(i)
        ystar <- rep(0, length(Y))
        for (k in 1:nind) { #### on retrace les effets al?atoires
          indiv <- which(id == unique(id)[k])
          ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }
        tree <- rpart::rpart(ystar ~ ., X) ### on construit l'arbre

        fhat <- predict(tree, X) #### pr?diction avec l'arbre
        for (k in 1:nind) { ### calcul des effets al?atoires par individu
          indiv <- which(id == unique(id)[k])
          V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv))
          btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat, id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde <- bay(btilde, Btilde, Z, id, sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y, fhat, Z[, , drop = FALSE], time, id, Btilde, 0, sigmahat, "none"))
        if (i > 1) inc <- (Vrai[i - 1] - Vrai[i]) / Vrai[i - 1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, Vraisemblance = Vrai, id = id, time = time)
          class(sortie) <- "APRF"
          return(sortie)
        }
      }
      sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, Vraisemblance = Vrai, id = id, time = time)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }

  for (i in 1:iter) {
    print(i)
    ystar <- rep(0, length(Y))
    for (k in 1:nind) { #### on retrace les effets al?atoires
      indiv <- which(id == unique(id)[k])
      ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }

    tree <- rpart::rpart(ystar ~ ., as.data.frame(X)) ### on construit l'arbre
    fhat <- predict(tree, as.data.frame(X)) #### pr?diction avec l'arbre
    for (k in 1:nind) { ### calcul des effets al?atoires par individu
      indiv <- which(id == unique(id)[k])
      K <- sto_analysis(sto, time[indiv])
      V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv)) + sigma2 * K
      btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
      omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ])
      epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }
    #### pr?diction du processus stochastique:
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat, id, Z, epsilonhat, Btilde, time, sigma2, sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde <- bay_sto(btilde, Btilde, Z, id, sigm, time, sigma2, sto) #### MAJ des param?tres de la variance des effets al?atoires.
    ### MAJ de la volatilit? du processus stochastique
    sigma2 <- gam_sto(sigm, id, Z, B, time, sigma2, sto, omega)
    Vrai <- c(Vrai, logV(Y, fhat, Z[, , drop = FALSE], time, id, Btilde, sigma2, sigmahat, sto))
    if (i > 1) inc <- (Vrai[i - 1] - Vrai[i]) / Vrai[i - 1]
    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_omega = id_omega, id_btilde = unique(id), omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }
  sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), id_omega = id_omega, omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id)
  class(sortie) <- "APRF"
  return(sortie)
}
#' Title
#'
#'
#' @import stats
#' @importFrom rpart predict
#'
#' @keywords internal
REEMtree <- function(X, Y, id, Z, iter = 10, time, sto, delta = 0.001) {
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0, nind, q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1, q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0, length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  id_omega <- sto

  if (class(sto) == "character") {
    if (sto == "none") {
      for (i in 1:iter) {
        print(i)
        ystar <- rep(0, length(Y))
        for (k in 1:nind) { #### on retrace les effets al?atoires
          indiv <- which(id == unique(id)[k])
          ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }

        tree <- rpart::rpart(ystar ~ ., X,control=rpart::rpart.control(cp=0)) ### on construit l'arbre
        Phi <- matrix(0, length(Y), length(unique(tree$where)))
        feuilles <- predict(tree, X)
        leaf <- unique(feuilles)
        for (p in 1:length(leaf)) {
          w <- which(feuilles == leaf[p])
          Phi[unique(w), p] <- 1
        }

        beta <- Moy(id, Btilde, sigmahat, Phi, Y, Z) ### fit des feuilles

        for (k in 1:length(unique(tree$where))) {
          ou <- which(tree$frame[, 5] == leaf[k])
          lee <- which(tree$frame[, 1] == "<leaf>")
          w <- intersect(ou, lee)
          tree$frame[w, 5] <- beta[k]
        }
        fhat <- predict(tree, as.data.frame(X))
        for (k in 1:nind) { ### calcul des effets al?atoires par individu
          indiv <- which(id == unique(id)[k])
          V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv))
          btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat, id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde <- bay(btilde, Btilde, Z, id, sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y, fhat, Z[, , drop = FALSE], time, id, Btilde, 0, sigmahat, sto))
        if (i > 1) inc <- (Vrai[i - 1] - Vrai[i]) / Vrai[i - 1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, vraisemblance = Vrai, id = id, time = time)
          class(sortie) <- "APRF"
          return(sortie)
        }
      }
      sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), sto = sto, Vraisemblance = Vrai, time = time, id = id)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }
  for (i in 1:iter) {
    print(i)
    ystar <- rep(0, length(Y))
    for (k in 1:nind) { #### on retrace les effets al?atoires
      indiv <- which(id == unique(id)[k])
      ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }

    tree <- rpart::rpart(ystar ~ ., X)
    Phi <- matrix(0, length(Y), length(unique(tree$where)))
    feuilles <- predict(tree, X)
    leaf <- unique(feuilles)
    for (p in 1:length(leaf)) {
      w <- which(feuilles == leaf[p])
      Phi[unique(w), p] <- 1
    }

    beta <- Moy_sto(id, Btilde, sigmahat, Phi, Y, Z, sto, time, sigma2) ### fit des feuilles

    for (k in 1:length(unique(tree$where))) {
      ou <- which(tree$frame[, 5] == leaf[k])
      lee <- which(tree$frame[, 1] == "<leaf>")
      w <- intersect(ou, lee)
      tree$frame[w, 5] <- beta[k]
    }
    fhat <- predict(tree, X)
    for (k in 1:nind) { ### calcul des effets al?atoires par individu
      indiv <- which(id == unique(id)[k])
      K <- sto_analysis(sto, time[indiv])
      V <- Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv)) + sigma2 * K
      btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
      omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k, ] - omega[indiv]
    }
    #### pr?diction du processus stochastique:
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat, id, Z, epsilonhat, Btilde, time, sigma2, sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde <- bay_sto(btilde, Btilde, Z, id, sigm, time, sigma2, sto) #### MAJ des param?tres de la variance des effets al?atoires.
    sigma2 <- gam_sto(sigm, id, Z, B, time, sigma2, sto, omega)
    Vrai <- c(Vrai, logV(Y, fhat, Z[, , drop = FALSE], time, id, Btilde, sigma2, sigmahat, sto))
    if (i > 1) inc <- (Vrai[i - 1] - Vrai[i]) / Vrai[i - 1]
    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), id_omega = id_omega, omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id)
      class(sortie) <- "APRF"
      return(sortie)
    }
  }
  sortie <- list(forest = tree, random_effects = btilde, var_random_effects = Btilde, sigma = sigmahat, id_btilde = unique(id), id_omega = id_omega, omega = omega, sigma_sto = sigma2, time = time, sto = sto, Vraisemblance = Vrai, id = id)
  class(sortie) <- "APRF"
  return(sortie)
}
#' Title
#'
#'
#'
#' @import stats
#'
#' @keywords internal
logV <- function(Y, f, Z, time, id, B, gamma, sigma, sto) {
  Vraisem <- 0
  if (sto == "none") {
    for (i in 1:length(unique(id))) {
      w <- which(id == unique(id)[i])
      V <- Z[w, , drop = FALSE] %*% B %*% t(Z[w, , drop = FALSE]) + diag(as.numeric(sigma), length(w), length(w))
      Vraisem <- Vraisem + log(det(V)) + t(Y[w] - f[w]) %*% solve(V) %*% (Y[w] - f[w])
    }
    return(Vraisem)
  }
  for (i in 1:length(unique(id))) {
    w <- which(id == unique(id)[i])
    K <- sto_analysis(sto, time[w])
    V <- Z[w, , drop = FALSE] %*% B %*% t(Z[w, , drop = FALSE]) + gamma * K + diag(as.numeric(sigma), length(w), length(w))
    Vraisem <- Vraisem + log(det(V)) + t(Y[w] - f[w]) %*% solve(V) %*% (Y[w] - f[w])
  }
  return(Vraisem)
}


