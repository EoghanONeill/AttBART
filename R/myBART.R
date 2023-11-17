#' @import Rcpp
#' @import progress
#' @import collapse
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom rmutil 'ddoublepois'
#' @useDynLib AttBART, .registration = TRUE
#' @export

attBart_no_w <- function(Xtrain,
                         ytrain,
                         m = 5,
                         node_min_size = 5, # Needs to be at least 1!
                         alpha = 0.95,
                         beta = 2,
                         nu = 3,
                         sigquant = 0.90,
                         k = 2,
                         lambda = NA,
                         sigest = NA,
                         sigmaf = NA,
                         n_burn = 1000,
                         n_post = 1000, # Number of observations post burn-in and thinning
                         n_thin = 1,
                         trans_prob = c(2.5, 2.5, 4) / 9, # Probabilities to grow, prune or change, respectively
                         max_bad_trees = 10,
                         sparse = FALSE, # The feature weighting of the DART model
                         a = 0.5, # ????
                         b = 1, # ????
                         seed = NA,
                         feature_weighting = FALSE,
                         sq_num_features = TRUE,
                         sq_ydiff_sigmu = TRUE,
                         centre_y = TRUE,
                         const_tree_weights = FALSE,
                         splitprob_as_weights = FALSE) { # Simple feature weighting
  if (!is.na(seed)) set.seed(seed)


  if(sparse == FALSE){
    splitprob_as_weights <- FALSE
  }

  if(splitprob_as_weights & feature_weighting){
    stop(" Cannot have both feature_weighting and splitprob_as_weights equal to TRUE")

  }

  # Transform y and X
  X_scaled <- scale(Xtrain)
  X_center <- attr(X_scaled, "scaled:center")
  X_scale <- attr(X_scaled, "scaled:scale")


  y_sd <- sd(ytrain)
  y_mean <- mean(ytrain)

  y_scale = (ytrain - y_mean)/y_sd

  # y_min <- min(y_scale)
  # y_max <- max(y_scale)

  if(centre_y){
    y_max <- max(y_scale)
    y_min <- min(y_scale)
  }else{
    y_max <- 0
    y_min <- 0
  }

  # Other variables
  sigma2 <- 1                          # !!!!!!!!!!!!!!
  mu_mu <- 0 # (y_min + y_max) / (2 * m)

  y_scale <- y_scale - (y_max + y_min)/2
  # tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))

  if(sq_ydiff_sigmu){
    # sigma2_mu <- ((max(y_scale)-min(y_scale))/(2 * k * sqrt(m)))^2
    sigma2_mu <- ((max(y_scale)-min(y_scale))*sqrt(m)/(2 * k))^2
  }else{
    # sigma2_mu <- (max(y_scale)-min(y_scale))/((2 * k * sqrt(m))^2)
    sigma2_mu <- (max(y_scale)-min(y_scale))*sqrt(m)/((2 * k)^2)
  }

  # sigma2_mu <- ((y_max - y_min) / (2 * k * sqrt(m)))^2

  alpha_s <- 1


  n <- length(y_scale)
  p <- ncol(X_scaled)
  s <- rep(1 / p, p) # probability vector to be used during the growing process for DART feature weighting
  rho <- p # For DART

  if (is.na(sigest)) {
    if (p < n) {
      df <- data.frame(X_scaled, y_scale)
      lmf <- lm(y_scale ~ ., df)
      sigest <- summary(lmf)$sigma
    } else {
      sigest <- sd(y_scale)
    }
  }
  qchi <- qchisq(1.0 - sigquant, nu)
  lambda <- (sigest * sigest * qchi) / nu # lambda parameter for sigma prior

  # if (is.na(sigmaf)) {
  #   tau <- (max(y_scale) - min(y_scale)) / (2 * k * sqrt(m))
  # } else {
  #   tau <- sigmaf / sqrt(m)
  # }
  sigma2 <- sigest^2
  tau <- 1

  # Total number of MCMC iterations
  n_iter <- n_burn + n_post * n_thin
  runif_matrix <- matrix(runif(m * n_iter), nrow = n_iter, ncol = m)

  # Storage containers
  tree_store <- vector(mode = "list", length = n_post)
  sigma2_store <- rep(NA, n_post)
  y_hat_store <- matrix(NA, ncol = n, nrow = n_post)
  att_weights_store <- vector(mode = "list", length = n_post)
  var_count <- rep(0, p)
  var_count_store <- matrix(0, ncol = p, nrow = n_post)
  s_prob_store <- matrix(0, ncol = p, nrow = n_post)
  alpha_MH_store <- matrix(NA, ncol = m, nrow = n_post)
  type_store <- matrix(NA, ncol = m, nrow = n_post)
  chosen_type_store <- matrix(NA, ncol = m, nrow = n_post)

  # Initialise trees using stumps
  curr_trees <- create_stumps(
    m = m,
    y = y_scale,
    X = X_scaled
  )

  # Set up progress bar
  # pb <- progress_bar$new(total = n_iter, format = "MCMC iterations [:bar] :current/:total in :elapsedfull, ETA: :eta")
  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = n_iter,
                             style = 3, width = 60,
                             title = 'Running attBART...')

  # MCMC iterations loop
  for (i in 1:n_iter) {
    utils::setTxtProgressBar(pb, i)

    # Loop through trees
    for (j in 1:m) {
      # We need the new and old trees for the likelihoods
      new_trees <- curr_trees

      # # Obtain the type
      # if (nrow(curr_trees[[j]]$tree_matrix) == 1) {
      #   type <- "grow"
      # } else {
      #   type <- sample(c("grow", "prune", "change"), 1, prob = trans_prob)
      # }


      type = sample_move(curr_trees[[j]], i, 100, #n_burn
                         trans_prob)

      # Generate a new tree based on the current
      new_trees[[j]] <- update_tree(
        y = y_scale,
        X = X_scaled,
        type = type,
        curr_tree = curr_trees[[j]],
        node_min_size = node_min_size,
        s = s,
        max_bad_trees = max_bad_trees
      )

      # # Feature weighting
      # if (feature_weighting) {
      #   tree_splitvars <- curr_trees[[j]]$tree_matrix[, "split_variable"]
      #   if (any(!is.na(tree_splitvars))) {
      #     # Remove the NAs
      #     splitvars <- sort(as.numeric(na.omit(tree_splitvars)))
      #   } else {
      #     # Else, return all variables
      #     splitvars <- 1:p
      #   }
      # } else {
      #   splitvars <- NA
      # }

      # Performing the Gibbs sampler:

      # 1. MH step --------------------------------------------------------------

      # Calculate the attention weights and the log likelihood using both the current trees and the proposed trees
      # (a) Calculations using the current trees

      if(const_tree_weights){
        att_weights_current <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
      }else{
        att_weights_current <- get_attention_no_w(curr_trees, X_scaled, tau, feature_weighting, sq_num_features,
                                                  splitprob_as_weights, s)
      }

      if(any(is.na(att_weights_current))){
        print("att_weights_current = ")
        print(att_weights_current)
        stop("attention weights contain nA")
      }

      if(any(  abs(1 - rowSums(att_weights_current)) > 0.0001     )){
        print("att_weights_current = ")
        print(att_weights_current)

        print("rowSums(att_weights_current) = ")
        print(rowSums(att_weights_current))

        print("which(rowSums(att_weights_current) != 1) = ")
        print(which(rowSums(att_weights_current) != 1))


        print("1 - rowSums(att_weights_current) = ")
        print(1 - rowSums(att_weights_current))

        stop("attention weights do not sum to 1")
      }


      # Create partial residuals conditional on the current tree
      no_j <- c(1:m)[-j]
      curr_partial_resid <- y_scale
      for (tree_ind in no_j) {
        curr_partial_resid <- curr_partial_resid - att_weights_current[, tree_ind] * get_prediction_no_w(curr_trees[[tree_ind]], X_scaled)
      }
      curr_partial_resid_rescaled <- curr_partial_resid / att_weights_current[, j]


      # (b) Calculations using the proposed tree
      if(const_tree_weights){
        att_weights_new <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
      }else{
        att_weights_new <- get_attention_no_w(new_trees, X_scaled, tau, feature_weighting, sq_num_features,
                                              splitprob_as_weights, s)
      }

      if(any(is.na(att_weights_new))){
        print("att_weights_new = ")
        print(att_weights_new)

        print("s = ")
        print(s)


        stop("attention weights contain nA")
      }

      if(any(  abs(1 - rowSums(att_weights_current)) > 0.0001  )){
        print("att_weights_new = ")
        print(att_weights_new)

        print("rowSums(att_weights_new) = ")
        print(rowSums(att_weights_new))

        print("which(rowSums(att_weights_new) != 1) = ")
        print(which(rowSums(att_weights_new) != 1))


        print("1 - rowSums(att_weights_new) = ")
        print(1 - rowSums(att_weights_new))


        stop("attention weights do not sum to 1")
      }
      # Create partial residuals conditional on the proposed tree
      no_j <- c(1:m)[-j]
      new_partial_resid <- y_scale
      for (tree_ind in no_j) {
        new_partial_resid <- new_partial_resid - att_weights_new[, tree_ind] * get_prediction_no_w(new_trees[[tree_ind]], X_scaled)
      }
      new_partial_resid_rescaled <- new_partial_resid / att_weights_new[, j]


      # (c) Obtain the Metropolis-Hastings probability
      curr_tree <- curr_trees[[j]]
      new_tree <- new_trees[[j]]
      alpha_MH <- get_MH_probability(
        X = X_scaled, curr_tree, new_tree,
        att_weights_current[, j], att_weights_new[, j],
        curr_partial_resid_rescaled, new_partial_resid_rescaled,
        type, trans_prob,
        alpha, beta,
        mu_mu, sigma2_mu, sigma2,
        node_min_size
      )

      # Accept new tree with probability alpha. Save additional information
      chosen_type <- type
      # if ((runif_matrix[i, j] < alpha_MH) | i < 4 ) {
      if ( runif_matrix[i, j] < alpha_MH  ) {
        curr_trees[[j]] <- new_tree
        curr_partial_resid_rescaled <- new_partial_resid_rescaled
        att_weights_current <- att_weights_new

        if (type == "grow") {
          var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] + 1
        } else if (type == "prune") {
          var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] - 1
        } else {
          var_count[curr_trees[[j]]$var[1]] <- var_count[curr_trees[[j]]$var[1]] - 1 # What if change step returned $var equal to c(0,0) ????
          var_count[curr_trees[[j]]$var[2]] <- var_count[curr_trees[[j]]$var[2]] + 1
        }
      } else {
        type <- "reject"
      }

      # Store the alphas if at the right place
      if ((i > n_burn) & ((i - n_burn) %% n_thin) == 0) {
        curr <- (i - n_burn) / n_thin
        alpha_MH_store[curr, j] <- alpha_MH
        type_store[curr, j] <- type
        chosen_type_store[curr, j] <- chosen_type
      }

      # 2. Update/Sample mu -----------------------------------------------------
      curr_trees[[j]] <- att_simulate_mu(curr_trees[[j]], curr_partial_resid_rescaled, att_weights_current[, j], mu_mu, sigma2_mu, sigma2)
    } # End loop through trees

    # 3. Update/Sample sigma2 ---------------------------------------------------
    y_hat <- rep(0, m)
    for (j in 1:m) {
      y_hat <- y_hat + get_prediction_no_w(curr_trees[[j]], X_scaled) * att_weights_current[, j]
    }
    sum_of_squares <- sum((y_scale - y_hat)^2)

    sigma2 <- update_sigma2(S = sum_of_squares, n, nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor q in 1:p is used to create new terminal nodes
    if (sparse & i > floor(n_iter * 0.1)) {
      s <- update_s(var_count, p, alpha_s)
    }

    # If at the right place, store everything
    if ((i > n_burn) & ((i - n_burn) %% n_thin) == 0) {
      curr <- (i - n_burn) / n_thin

      tree_store[[curr]] <- curr_trees
      sigma2_store[curr] <- sigma2
      y_hat_store[curr, ] <- y_hat
      var_count_store[curr, ] <- var_count
      s_prob_store[curr, ] <- s
      att_weights_store[[curr]] <- att_weights_current
    }
    # pb$tick()
  } # End loop through MCMC iterations

  cat("\n") # Make sure progress bar ends on a new line

  return(list(
    trees = tree_store,
    sigma2 = sigma2_store * y_sd^2,
    y_hat =  (y_hat_store+(y_max + y_min)/2)*y_sd + y_mean,
    var_count_store = var_count_store,
    s = s_prob_store,
    center = X_center,
    scale = X_scale,
    scaledtrainingdata = X_scaled,
    MH_prob = alpha_MH_store,
    att_weights = att_weights_store,
    types = type_store,
    chosen_types = chosen_type_store,
    npost = n_post,
    nburn = n_burn,
    nthin = n_thin,
    ntrees = m,
    y_sd = y_sd,
    y_mean = y_mean,
    feature_weighting = feature_weighting,
    tau = tau,
    y_max = y_max,
    y_min = y_min,
    const_tree_weights = const_tree_weights,
    sq_num_features = sq_num_features,
    splitprob_as_weights = splitprob_as_weights
  ))
}
