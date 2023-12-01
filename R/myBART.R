



make_01_norm <- function(x) {
  a <- min(x)
  b <- max(x)
  return(function(y0) (y0 - a) / (b - a))
}

make_normalized <- function(x) {
  a <- mean(x)
  b <- sd(x)
  if(b == 0 | is.na(b)){
    b <- 1
  }
  return(function(y0) (y0 - a) / (b))
}

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
                         alpha_a = 0.5, # Linero alpha prior parameter
                         alpha_b = 1, # Linero alpha prior parameter
                         seed = NA,
                         feature_weighting = FALSE,
                         sq_num_features = TRUE,
                         sq_ydiff_sigmu = TRUE,
                         centre_y = TRUE,
                         const_tree_weights = FALSE,
                         splitprob_as_weights = FALSE,
                         tau = 1,
                         alpha_prior = FALSE,
                         sigma_mu_prior = FALSE,
                         update_tau = FALSE,
                         covariate_scaling = "normalize",
                         warm_weight_start = 0,
                         include_w = FALSE,
                         fix_epsilon_w = TRUE,
                         epsilon_w = 0.5,
                         w_prior = "normal",
                         tau_w_hyperprior = FALSE) { # Simple feature weighting


  if(include_w & !fix_epsilon_w ){
    stop("Currently code is only written for fixed epsilon")
  }

  if(include_w & const_tree_weights ){
    stop("Cannot include w and have constant tree weights. Maybe constant attention weights and varying w can be implemented in future implementations.")
  }


  if(!(w_prior %in% c("normal", "Dirichlet"))){
    stop("w_prior must be 'normal' or 'Dirichlet'.")
  }

  if (!is.na(seed)) set.seed(seed)

  if(!(covariate_scaling %in% c("none", "ECDF", "normalize"))){
    stop("covariate_scaling must be none, ECDF, or normalize")
  }

  if(sparse == FALSE){
    splitprob_as_weights <- FALSE
  }

  if(splitprob_as_weights & feature_weighting){
    stop(" Cannot have both feature_weighting and splitprob_as_weights equal to TRUE")

  }

  # # Transform y and X
  # X_scaled <- scale(Xtrain)
  # X_center <- attr(X_scaled, "scaled:center")
  # X_scale <- attr(X_scaled, "scaled:scale")

  Xtrain <- as.matrix(Xtrain)

  X_scaled <- matrix(NA,
                     nrow = nrow(Xtrain),
                     ncol = ncol(Xtrain))

  scale_x_funcs   <- list()
  for(i in 1:ncol(Xtrain)) {

    if(covariate_scaling == "none") scale_x_funcs[[i]] <- identity
    if(covariate_scaling == "ECDF") scale_x_funcs[[i]] <- ecdf(Xtrain[,i])
    if(covariate_scaling == "normalize") scale_x_funcs[[i]] <- make_normalized(Xtrain[,i])

    if(length(unique(Xtrain[,i])) == 1) scale_x_funcs[[i]] <- identity
    if(length(unique(Xtrain[,i])) == 2) scale_x_funcs[[i]] <- make_01_norm(Xtrain[,i])
  }
  for(i in 1:ncol(Xtrain)) {
    X_scaled[,i] <- scale_x_funcs[[i]](Xtrain[,i])
  }



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
    sigma2_mu <- (max(y_scale)-min(y_scale))*m/((2 * k)^2)
  }

  # sigma2_mu <- ((y_max - y_min) / (2 * k * sqrt(m)))^2



  n <- length(y_scale)
  p <- ncol(X_scaled)
  s <- rep(1 / p, p) # probability vector to be used during the growing process for DART feature weighting
  rho <- p # For DART

  alpha_s <- 1 # p

  alpha_scale <- p



  tau_rate <- 10

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
  # tau <- 1

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

  tau_store <- rep(NA, n_post)


  # Initialise trees using stumps
  curr_trees <- create_stumps(
    m = m,
    y = y_scale,
    X = X_scaled
  )

  # if(const_tree_weights){
  #   att_weights_current <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
  # }else{
  #   att_weights_current <- get_attention_no_w(curr_trees, X_scaled, tau, feature_weighting, sq_num_features,
  #                                             splitprob_as_weights, s)
  # }








  if(const_tree_weights | (warm_weight_start > 0)
     ){
    # att_weights_current <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)

    att_weights_current_unnorm <- matrix( 1, nrow = nrow(X_scaled), ncol = m)
    att_weights_current_denoms <- rowSums(att_weights_current_unnorm)
    att_weights_current <- att_weights_current_unnorm/att_weights_current_denoms

    # print("initial weights unorm = ")
    # print(att_weights_current_unnorm)
    # print("initial weights = ")
    # print(att_weights_current)

    att_weights_new_unnorm <- att_weights_current_unnorm
    att_weights_new_denoms <- att_weights_current_denoms
    att_weights_new <- att_weights_current

  }else{
    att_weights_current_unnorm <- get_unnorm_att_all_no_w(curr_trees, X_scaled, tau, feature_weighting, sq_num_features,
                                              splitprob_as_weights, s, FALSE)
    att_weights_current_denoms <- rowSums(att_weights_current_unnorm)
    att_weights_current <- att_weights_current_unnorm/att_weights_current_denoms
#
#     print("initial weights unorm = ")
#     print(att_weights_current_unnorm)
#     print("initial weights = ")
#     print(att_weights_current)

    att_weights_new_unnorm <- att_weights_current_unnorm
    att_weights_new_denoms <- att_weights_current_denoms
    att_weights_new <- att_weights_current
  }



  # the include_w option will be included separately to the ocnstant weights option for hte purpose of testing






  # Set up progress bar
  # pb <- progress_bar$new(total = n_iter, format = "MCMC iterations [:bar] :current/:total in :elapsedfull, ETA: :eta")
  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = n_iter,
                             style = 3, width = 60,
                             title = 'Running attBART...')



  treepredmat <- matrix(NA, nrow = n, ncol = m)


  if(warm_weight_start >= n_iter){
    stop("Number of warm start iterations is greater than or wqual to total number of MCMC iterations.")
  }


  for(j in 1:m){

    treepredmat[,j] <- get_prediction_no_w(curr_trees[[j]], X_scaled)
  }


  if(include_w){

    w_vec <- rep(1/m, m)
    w_prior_mean <- rep(1/m, m)
    tau_w <- 1

    a_w <- 1
    b_w <- 1


    if(w_prior == "Dirichlet"){
      xi_w_vec <- rep(1/m, m)
      alpha_w_par <- m

      a_w_vec <- alpha_w_par*xi_w_vec
      v_w_vec <- rep(0,m)

      a_step_temp <- -(m/alpha_w_par)
      m_step_temp <- 0.5*a_step_temp
      s_step_temp <- -m_step_temp*tan((pi/10)/(m/alpha_w_par))

      # print("a_step_temp = ")
      # print(a_step_temp)
      #
      # print("m_step_temp = ")
      # print(m_step_temp)
      #
      # print("s_step_temp = ")
      # print(s_step_temp)

    }

    Adjusted_att_weights_current <- (att_weights_current*(1 - epsilon_w) ) %r+% (w_vec*epsilon_w)
    Adjusted_att_weights_new <- (att_weights_new*(1 - epsilon_w)) %r+% (w_vec*epsilon_w)

    # Adjusted_att_weights_current <- Adjusted_att_weights_current/rowSums(Adjusted_att_weights_current)
    # Adjusted_att_weights_new <- Adjusted_att_weights_new/rowSums(Adjusted_att_weights_new)

    w_by_pred_sums <- treepredmat %*% w_vec #rowSums(treepredmat %r+% w_vec )


    wvec_store <- matrix(NA, ncol = m, nrow = n_post)
    tau_w_store <- rep(NA,n_post)


  }

  # print("treepredmat[1,] %*% w_vec")
  # print(treepredmat[1,] %*% w_vec)
  #
  # print("treepredmat[1,] %*% w_vec")
  # print(treepredmat[1,] %*% w_vec)

  y_hat_unnorm <- rowSums(treepredmat*att_weights_current_unnorm)


  # print('cbind(w_by_pred_sums,y_hat_unnorm/att_weights_current_denoms, y_scale )')
  # print(cbind(w_by_pred_sums,y_hat_unnorm/att_weights_current_denoms, y_scale ))


  ######### MCMC iterations loop ##########
  for (i in 1:n_iter) {
    utils::setTxtProgressBar(pb, i)

    if(i == warm_weight_start +1){
      att_weights_current_unnorm <- get_unnorm_att_all_no_w(curr_trees, X_scaled, tau, feature_weighting, sq_num_features,
                                                            splitprob_as_weights, s, FALSE)
      att_weights_current_denoms <- rowSums(att_weights_current_unnorm)
      att_weights_current <- att_weights_current_unnorm/att_weights_current_denoms

      att_weights_new_unnorm <- att_weights_current_unnorm
      att_weights_new_denoms <- att_weights_current_denoms
      att_weights_new <- att_weights_current
    }

    if(include_w){


      Adjusted_att_weights_current <- (att_weights_current*(1 - epsilon_w)) %r+% (w_vec*epsilon_w)
      Adjusted_att_weights_new <- (att_weights_new*(1 - epsilon_w)) %r+% (w_vec*epsilon_w)

      # Adjusted_att_weights_current <- Adjusted_att_weights_current/rowSums(Adjusted_att_weights_current)
      # Adjusted_att_weights_new <- Adjusted_att_weights_new/rowSums(Adjusted_att_weights_new)

    }

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



      type = sample_move(curr_trees[[j]], i, 0, #n_burn
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

      # if(const_tree_weights){
      #   att_weights_current <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
      # }else{
      #   att_weights_current <- get_attention_no_w(curr_trees, X_scaled, tau, feature_weighting, sq_num_features,
      #                                             splitprob_as_weights, s)
      # }
      #
      #
      # att_vec_current <- get_unnorm_att_1tree_no_w(curr_trees[[j]], X_scaled, tau, feature_weighting, sq_num_features,
      #                           splitprob_as_weights, s)


      # if(any(is.na(att_weights_current))){
      #   print("att_weights_current = ")
      #   print(att_weights_current)
      #   stop("attention weights contain nA")
      # }
      #
      # if(any(  abs(1 - rowSums(att_weights_current)) > 0.0001     )){
      #   print("att_weights_current = ")
      #   print(att_weights_current)
      #
      #   print("rowSums(att_weights_current) = ")
      #   print(rowSums(att_weights_current))
      #
      #   print("which(rowSums(att_weights_current) != 1) = ")
      #   print(which(rowSums(att_weights_current) != 1))
      #
      #
      #   print("1 - rowSums(att_weights_current) = ")
      #   print(1 - rowSums(att_weights_current))
      #
      #   print("iteration i = ")
      #   print(i)
      #
      #   print("tree j = ")
      #   print(j)
      #
      #
      #   stop("attention weights do not sum to 1")
      # }


      # Create partial residuals conditional on the current tree
      # no_j <- c(1:m)[-j]
      # curr_partial_resid <- y_scale
      # for (tree_ind in no_j) {
      #   curr_partial_resid <- curr_partial_resid - att_weights_current[, tree_ind] * treepredmat[,tree_ind]#get_prediction_no_w(curr_trees[[tree_ind]], X_scaled)
      # }



      if(include_w){
        w_by_pred_sums_less_j <- w_by_pred_sums - w_vec[j] * treepredmat[,j]
        y_hat_unnorm_less_j <- (y_hat_unnorm - att_weights_current_unnorm[, j] * treepredmat[,j])

        curr_partial_resid <- y_scale - ((1-epsilon_w)*y_hat_unnorm_less_j /att_weights_current_denoms) - epsilon_w*w_by_pred_sums_less_j
        # new_partial_resid_rescaled <- new_partial_resid / Adjusted_att_weights_new[, j]

        Adjusted_att_weights_current_j <- att_weights_current[,j]*(1 - epsilon_w) + w_vec[j]*epsilon_w
        curr_partial_resid_rescaled <- curr_partial_resid /  Adjusted_att_weights_current_j

        # print("att_weights_current[1:3,]")
        # print(att_weights_current[1:3,])
        #
        # print("Adjusted_att_weights_current_j[1:10]")
        # print(Adjusted_att_weights_current_j[1:10])
        #
        #
        # print("curr_partial_resid[1:10]")
        # print(curr_partial_resid[1:10])
        #
        # print("curr_partial_resid_rescaled[1:10]")
        # print(curr_partial_resid_rescaled[1:10])
        #
        # print("mean(curr_partial_resid)")
        # print(mean(curr_partial_resid))
        #
        # print("mean(curr_partial_resid_rescaled)")
        # print(mean(curr_partial_resid_rescaled))
        #
        # print("sd(curr_partial_resid)")
        # print(sd(curr_partial_resid))
        #
        # print("sd(curr_partial_resid_rescaled)")
        # print(sd(curr_partial_resid_rescaled))


      }else{
        # curr_partial_resid <- y_scale - rowSums(att_weights_current[, -j] * treepredmat[,-j])
        curr_partial_resid <- y_scale - y_hat_unnorm/att_weights_current_denoms + att_weights_current[, j] * treepredmat[,j]

        # y_hat_unnorm_less_j <- (y_hat_unnorm - att_weights_current_unnorm[, j] * treepredmat[,j])
        # curr_partial_resid <- y_scale - y_hat_unnorm_less_j /att_weights_current_denoms

        curr_partial_resid_rescaled <- curr_partial_resid / att_weights_current[, j]

        # print("mean(curr_partial_resid)")
        # print(mean(curr_partial_resid))
        #
        # print("mean(curr_partial_resid_rescaled)")
        # print(mean(curr_partial_resid_rescaled))

      }

      # (b) Calculations using the proposed tree
      # if(const_tree_weights){
      #   att_weights_new <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
      # }else{
      #   att_weights_new <- get_attention_no_w(new_trees, X_scaled, tau, feature_weighting, sq_num_features,
      #                                         splitprob_as_weights, s)
      # }

      # att_weights_new_unnorm <- att_weights_current_unnorm

      if(const_tree_weights | (i <= warm_weight_start) ){
        # att_weights_new <- matrix(1/m, nrow = nrow(X_scaled), ncol = m)
      }else{
        att_weights_new_unnorm <- att_weights_current_unnorm

        att_weights_new_unnorm[,j] <- get_unnorm_att_1tree_no_w(new_trees[[j]], X_scaled, tau, feature_weighting, sq_num_features,
                                                              splitprob_as_weights, s, FALSE)
        att_weights_new_denoms <- att_weights_current_denoms - att_weights_current_unnorm[,j] + att_weights_new_unnorm[,j]
        att_weights_new <- att_weights_new_unnorm/att_weights_new_denoms

        if(any(abs(rowSums(att_weights_new) - 1)>0.01)){

          print("rowSums(att_weights_current) = ")
          print(rowSums(att_weights_current))
          print("att_weights_current = ")
          print(att_weights_current)

          print("rowSums(att_weights_new) = ")
          print(rowSums(att_weights_new))
          print("att_weights_new = ")
          print(att_weights_new)

          stop("weights do not sum to 1")
        }


      }





      # if(any(is.na(att_weights_new))){
      #   print("att_weights_new = ")
      #   print(att_weights_new)
      #
      #   print("s = ")
      #   print(s)
      #
      #
      #   stop("attention weights contain nA")
      # }

      # if(any(  abs(1 - rowSums(att_weights_current)) > 0.0001  )){
      #   print("att_weights_new = ")
      #   print(att_weights_new)
      #
      #   print("rowSums(att_weights_new) = ")
      #   print(rowSums(att_weights_new))
      #
      #   print("which(rowSums(att_weights_new) != 1) = ")
      #   print(which(rowSums(att_weights_new) != 1))
      #
      #
      #   print("1 - rowSums(att_weights_new) = ")
      #   print(1 - rowSums(att_weights_new))
      #
      #
      #   stop("attention weights do not sum to 1")
      # }
      # Create partial residuals conditional on the proposed tree
      # no_j <- c(1:m)[-j]
      # new_partial_resid <- y_scale
      # for (tree_ind in no_j) {
      #   new_partial_resid <- new_partial_resid - att_weights_new[, tree_ind] * treepredmat[,tree_ind]#get_prediction_no_w(new_trees[[tree_ind]], X_scaled)
      # }

      # new_partial_resid <- y_scale - rowSums(att_weights_new[, -j] * treepredmat[,-j])





      if(include_w){


        # Adjusted_att_weights_current <- att_weights_current*(1 - epsilon_w) %r+% w_vec*epsilon_w
        # Adjusted_att_weights_new <- att_weights_new*(1 - epsilon_w) %r+% w_vec*epsilon_w

        # Adjusted_att_weights_current <- Adjusted_att_weights_current/rowSums(Adjusted_att_weights_current)
        # Adjusted_att_weights_new <- Adjusted_att_weights_new/rowSums(Adjusted_att_weights_new)


        w_by_pred_sums_less_j <- w_by_pred_sums - w_vec[j] * treepredmat[,j]
        y_hat_unnorm_less_j <- (y_hat_unnorm - att_weights_current_unnorm[, j] * treepredmat[,j])

        new_partial_resid <- y_scale - ((1-epsilon_w)*y_hat_unnorm_less_j /att_weights_new_denoms) - epsilon_w*w_by_pred_sums_less_j
        # new_partial_resid_rescaled <- new_partial_resid / Adjusted_att_weights_new[, j]

        Adjusted_att_weights_new_j <- att_weights_new[,j]*(1 - epsilon_w) + w_vec[j]*epsilon_w
        new_partial_resid_rescaled <- new_partial_resid /  ( Adjusted_att_weights_new_j)


        # print("att_weights_new[1:3,]")
        # print(att_weights_new[1:3,])
        #
        # print("Adjusted_att_weights_new_j[1:10]")
        # print(Adjusted_att_weights_new_j[1:10])
        #
        #
        # print("new_partial_resid[1:10]")
        # print(new_partial_resid[1:10])
        #
        # print("new_partial_resid_rescaled[1:10]")
        # print(new_partial_resid_rescaled[1:10])
        #
        #
        # print("mean(new_partial_resid)")
        # print(mean(new_partial_resid))
        #
        # print("mean(new_partial_resid_rescaled)")
        # print(mean(new_partial_resid_rescaled))
        #
        # print("sd(new_partial_resid)")
        # print(sd(new_partial_resid))
        #
        # print("sd(new_partial_resid_rescaled)")
        # print(sd(new_partial_resid_rescaled))





      }else{

        y_hat_unnorm_less_j <- (y_hat_unnorm - att_weights_current_unnorm[, j] * treepredmat[,j])
        new_partial_resid <- y_scale - y_hat_unnorm_less_j /att_weights_new_denoms

        new_partial_resid_rescaled <- new_partial_resid / att_weights_new[, j]



      }

      # if(include_w){
      #   # there is probably a more efficient way of doing this
      #   Adjusted_att_weights_new <- (att_weights_new*(1 - epsilon_w)) %r+% (w_vec*epsilon_w)
      #   # Adjusted_att_weights_new <- Adjusted_att_weights_new/rowSums(Adjusted_att_weights_new)
      #
      # }

      # (c) Obtain the Metropolis-Hastings probability
      curr_tree <- curr_trees[[j]]
      new_tree <- new_trees[[j]]

      if((nrow(new_tree$tree_matrix) == nrow(curr_tree$tree_matrix) ) & (type != "change" )){
        alpha_MH <- 0
        # print("no good trees")
      }else{

        if(include_w){

          if(any(is.na(Adjusted_att_weights_new_j))){
            print("temp_curr_attweights = ")
            print(temp_curr_attweights)

            print("temp_new_attweights = ")
            print(temp_new_attweights)

            print("Adjusted_att_weights_current_j = ")
            print(Adjusted_att_weights_current_j)

            print("Adjusted_att_weights_new_j = ")
            print(Adjusted_att_weights_new_j)

            print("j =" )
            print(j)

            print("i =" )
            print(i)


            stop("NA attention weights")

          }

          temp_new_attweights <- abs(Adjusted_att_weights_new_j)
          temp_curr_attweights <- abs(Adjusted_att_weights_current_j)
        }else{
          temp_new_attweights <- att_weights_new[, j]
          temp_curr_attweights <- att_weights_current[, j]
        }

        alpha_MH <- get_MH_probability(
          X = X_scaled, curr_tree, new_tree,
          temp_curr_attweights, #att_weights_current[, j],
          temp_new_attweights, # att_weights_new[, j],
          curr_partial_resid_rescaled, new_partial_resid_rescaled,
          type, trans_prob,
          alpha, beta,
          mu_mu, sigma2_mu, sigma2,
          node_min_size
        )
        # print("calculated MH probability")
        # print("new_tree = ")
        # print(new_tree)
      }


      # print("alpha_MH = ")
      # print(alpha_MH)

      if(is.na(alpha_MH)){
        print("alpha_MH = ")
        print(alpha_MH)

        print("curr_partial_resid_rescaled = ")
        print(curr_partial_resid_rescaled)

        print("new_partial_resid_rescaled = ")
        print(new_partial_resid_rescaled)

        print("temp_curr_attweights = ")
        print(temp_curr_attweights)

        print("temp_new_attweights = ")
        print(temp_new_attweights)

        print("j =" )
        print(j)

        print("i =" )
        print(i)

        stop("NA alpha_MH")

      }



      # y_hat_unnorm <- y_hat_unnorm - att_weights_current_unnorm[, j] * treepredmat[,j]


      # Accept new tree with probability alpha. Save additional information
      chosen_type <- type
      # if ((runif_matrix[i, j] < alpha_MH) | i < 4 ) {
      if ( runif_matrix[i, j] < alpha_MH  ) {
        curr_trees[[j]] <- new_tree
        curr_partial_resid_rescaled <- new_partial_resid_rescaled
        att_weights_current <- att_weights_new
        if( (!const_tree_weights) & (i > warm_weight_start)){
          att_weights_current_unnorm <- att_weights_new_unnorm
          att_weights_current_denoms <- att_weights_new_denoms
        }
        if (type == "grow") {
          var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] + 1
        } else if (type == "prune") {
          var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] - 1
        } else {
          var_count[curr_trees[[j]]$var[1]] <- var_count[curr_trees[[j]]$var[1]] - 1 # What if change step returned $var equal to c(0,0) ????
          var_count[curr_trees[[j]]$var[2]] <- var_count[curr_trees[[j]]$var[2]] + 1
        }


        if(include_w){
          Adjusted_att_weights_current_j <- Adjusted_att_weights_new_j
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

      if(include_w){
        # temp_new_attweights <- abs(Adjusted_att_weights_new_j)
        temp_curr_attweights <- abs(Adjusted_att_weights_current_j)
      }else{
        # temp_new_attweights <- att_weights_new[, j]
        temp_curr_attweights <- att_weights_current[, j]
      }


      # curr_trees[[j]] <- att_simulate_mu(curr_trees[[j]], curr_partial_resid_rescaled, att_weights_current[, j], mu_mu, sigma2_mu, sigma2)
      curr_trees[[j]] <- att_simulate_mu(curr_trees[[j]], curr_partial_resid_rescaled, temp_curr_attweights, mu_mu, sigma2_mu, sigma2)
      treepredmat[,j] <- get_prediction_no_w(curr_trees[[j]], X_scaled)

      # y_hat_unnorm <- y_hat_unnorm_less_j + treepredmat[,j] * att_weights_current_unnorm[, j]
      y_hat_unnorm <- rowSums(treepredmat * att_weights_current_unnorm)

      if(include_w){
        # w_by_pred_sums <- w_by_pred_sums_less_j + treepredmat[,j] * w_vec[j]
        w_by_pred_sums <- treepredmat %*% w_vec  # rowSums(treepredmat %r*% w_vec ) # rowSums(treepredmat %r+% w_vec )

        # print('treepredmat[1:3,] = ')
        # print(treepredmat[1:3,])
        # print('w_by_pred_sums[1:5] = ')
        # print(w_by_pred_sums[1:5])


      }


    } # End loop through trees

    ############# 3. Update/Sample sigma2 ---------------------------------------------------
    # y_hat <- rep(0, n)
    # for (j in 1:m) {
    #   # y_hat <- y_hat + get_prediction_no_w(curr_trees[[j]], X_scaled) * att_weights_current[, j]
    #   y_hat <- y_hat + treepredmat[,j] * att_weights_current[, j]
    # }


    if(include_w){
      y_hat <-  (1 - epsilon_w) * (y_hat_unnorm/att_weights_current_denoms) + epsilon_w*w_by_pred_sums
    }else{
      y_hat <- y_hat_unnorm/att_weights_current_denoms   # rowSums(treepredmat*att_weights_current)
    }

    sum_of_squares <- sum((y_scale - y_hat)^2)

    # if(sum_of_squares > 1000){
    #   print('cbind(y_scale, y_hat) = ')
    #   print(cbind(y_scale, y_hat))
    #
    #   print('y_hat_unnorm/att_weights_current_denoms = ')
    #   print(y_hat_unnorm/att_weights_current_denoms)
    #
    #   print('w_by_pred_sums= ')
    #   print(w_by_pred_sums)
    #
    #   print("sigma2 = ")
    #   print(sigma2)
    #
    #   print("i =" )
    #   print(i)
    #
    #   stop("big sum_of_squares")
    #
    # }



    sigma2 <- update_sigma2(S = sum_of_squares, n, nu, lambda)

    # print("w_vec = ")
    # print(w_vec)
    #
    # print("sigma2 = ")
    # print(sigma2)


    if(sigma2 > 100){
      print('cbind(y_scale, y_hat)[1:10,] = ')
      print(cbind(y_scale, y_hat)[1:10,])

      print(' ((1- epsilon_w) * y_hat_unnorm/att_weights_current_denoms)[1:5] = ')
      print(((1- epsilon_w) * y_hat_unnorm/att_weights_current_denoms)[1:5])

      print('epsilon_w*w_by_pred_sums[1:5]= ')
      print(epsilon_w*w_by_pred_sums[1:5])

      print("w_vec = ")
      print(w_vec)



      print("treepredmat[1:10,] = ")
      print(treepredmat[1:10,])


      print("sigma2 = ")
      print(sigma2)

      print("j =" )
      print(j)

      print("i =" )
      print(i)



      stop("big sigma2")

    }
    ######### update s and alpha_s ##########################

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor q in 1:p is used to create new terminal nodes
    if (sparse & (i > floor(n_burn * 0.5))) {
      s <- update_s(var_count, p, alpha_s)
      if(alpha_prior){
        alpha_s <- update_alpha(s, alpha_scale, alpha_a, alpha_b)
      }
    }


    ####### update sigma_mu #########################
    if(sigma_mu_prior){
      sigma2_mu <-  update_sigma_mu(curr_trees, sigma2_mu)
    }



    ####### update tau ####################

    if( ((! const_tree_weights) &  (i > warm_weight_start))  & update_tau){
      prop_tau <- max(tau*(5^(runif(n = 1,min = -1,max = 1))), 0.0001)

      att_weights_new_unnorm <- get_unnorm_att_all_no_w(curr_trees, X_scaled, prop_tau, feature_weighting, sq_num_features,
                                                            splitprob_as_weights, s, FALSE)
      att_weights_new_denoms <- rowSums(att_weights_new_unnorm)
      att_weights_new <- att_weights_new_unnorm/att_weights_new_denoms

      y_hat_unnorm_new <- rep(0, n)
      for (j in 1:m) {
        # y_hat <- y_hat + get_prediction_no_w(curr_trees[[j]], X_scaled) * att_weights_current[, j]
        y_hat_unnorm_new <- y_hat_unnorm_new + treepredmat[,j] * att_weights_new_unnorm[, j]
      }
      y_hat_new <- y_hat_unnorm_new/att_weights_new_denoms
      # sum_of_squares_new <- sum((y_scale - y_hat_new)^2)

      lik_orig <- sum(dnorm(x = y_scale,
                                   mean = y_hat,
                                   sd = sqrt(sigma2), log = TRUE))

      lik_prop <- sum(dnorm(x = y_scale,
                        mean = y_hat_new,
                        sd = sqrt(sigma2), log = TRUE))


      # l_new <-   - sum_of_squares_new/(2*sigma2)  + dexp(prop_tau,tau_rate, log = TRUE) - log(tau)
      # l_old <-   -  sum_of_squares/(2*sigma2)  + dexp(tau,tau_rate, log = TRUE) - log(prop_tau)

      l_new <-   lik_prop  + dexp(prop_tau,tau_rate, log = TRUE) - log(tau)
      l_old <-   lik_orig  + dexp(tau,tau_rate, log = TRUE) - log(prop_tau)


      if(runif(1) < exp(l_new - l_old)){

        # print("tau accepted")
        # print(prop_tau)

        tau <- prop_tau

        att_weights_current_unnorm <- att_weights_new_unnorm
        att_weights_current_denoms <- att_weights_new_denoms
        att_weights_current <- att_weights_new

        y_hat <- y_hat_new
        y_hat_unnorm <- y_hat_unnorm_new


        # if(include_w){
        #   # maybe do not need to update adjusted weights here because will be updated after new draw of w
        # }

      }

    }


    ###### update w ################

    # add checks for number of draws and for warm start etc

    if(include_w  &  (i > warm_weight_start) ){

      if(w_prior == "normal"){
         # define scaled residuals
        resid_w <- (y_scale - (1 - epsilon_w) * y_hat_unnorm / att_weights_current_denoms)/epsilon_w


        # obtain parameters for full conditional draw
        # these calculations can probably be made more efficient

        # print("w_varmat = ")
        # print(w_varmat)

        # print("treepredmat = ")
        # print(treepredmat)
        # print("tau_w = ")
        # print(tau_w)
        # print(" epsilon_w^2 = ")
        # print( epsilon_w^2)
        # print("sigma2 = ")
        # print(sigma2)
        # print("(tau_w * epsilon_w^2)/sigma2 = ")
        # print((tau_w * epsilon_w^2)/sigma2)
        #
        #
        # print("diag(m)*(tau_w * epsilon_w^2)/sigma2 = ")
        # print(diag(m)*(tau_w * epsilon_w^2)/sigma2)

        w_varmat <- solve( t(treepredmat) %*% treepredmat + diag(m)*(tau_w * epsilon_w^2)/sigma2)

        # print("w_varmat = ")
        # print(w_varmat)
        #
        # print("treepredmat = ")
        # print(treepredmat)
        #
        # print("w_prior_mean = ")
        # print(w_prior_mean)


        w_mean <- w_varmat %*% ( t(treepredmat) %*% resid_w   + w_prior_mean*(tau_w * epsilon_w^2)/sigma2  )
        w_varmat <- sigma2*w_varmat

        # draw new w value
        w_vec <- as.vector(rmvnorm(n = 1, mean = w_mean , sigma = w_varmat))

        # update w_by_pred_sums and any other predictions that can be updated
        w_by_pred_sums <- treepredmat %*% w_vec # rowSums(treepredmat %r*% w_vec ) # rowSums(treepredmat %r+% w_vec )

        # draw new tau_w if it is an option

        if(tau_w_hyperprior){
          tau_w <- rgamma(n = 1,
                          shape = a_w + m/2,
                          rate = b_w +  sum((w_vec - w_prior_mean)^2)*epsilon_w^2 /(2*sigma2))
        }


        tau_w <- tau_w/sum(tau_w)

        if(any(is.na("w_vec")) | is.na("tau_w"))  {
          print("w_varmat = ")
          print(w_varmat)

          print("treepredmat = ")
          print(treepredmat)

          print("w_prior_mean = ")
          print(w_prior_mean)

          print("w_vec = ")
          print(w_vec)

          print("tau_w = ")
          print(tau_w)

          stop("NA in w vec or tau_w")

        }

      }


      if(w_prior == "Dirichlet"){


        # sample reparametrized v_j values



        for(j in 1:m){

          step_size <-  abs( a_step_temp*(0.5 - atan( ( v_w_vec[j] -  m_step_temp) / s_step_temp)/pi ) )

          # print("v_w_vec = ")
          # print(v_w_vec)
          #
          # print("step_size = ")
          # print(step_size)
          #
          # print("j = ")
          # print(j)

          v_j_prop <- rnorm(n = 1,
                            mean = v_w_vec[j],
                            sd = step_size
                            )

          l_prop_prob <- dnorm(x = v_j_prop,
                             mean = v_w_vec[j],
                             sd = step_size,
                             log = TRUE)


          reverse_step_size <-  abs(a_step_temp*(0.5 - atan( ( v_j_prop -  m_step_temp) / s_step_temp)/pi ))

          l_reverse_prob <- dnorm(x = v_w_vec[j],
                                mean = v_j_prop,
                                sd = reverse_step_size,
                                log = TRUE)


          # prior_prop <- exp(a_w_vec[j] * v_j_prop - exp(v_j_prop))/gamma(a_w_vec[j] )
          # prior_orig <- exp(a_w_vec[j] * v_w_vec[j] - exp(v_w_vec[j]))/gamma(a_w_vec[j] )

          l_prior_prop <- a_w_vec[j] * v_j_prop - exp(v_j_prop) - lgamma(a_w_vec[j] )
          l_prior_orig <- a_w_vec[j] * v_w_vec[j] - exp(v_w_vec[j]) -lgamma(a_w_vec[j] )

          prop_v_w_vec <- v_w_vec
          prop_v_w_vec[j] <- v_j_prop


          # print("v_j_prop = ")
          # print(v_j_prop)

          prop_w_vec <- exp(prop_v_w_vec)/sum(exp(prop_v_w_vec))

          # now require likelihoods from new and old w values
          w_by_pred_sums <- treepredmat %*% w_vec  # rowSums(treepredmat %r*% w_vec ) # rowSums(treepredmat %r+% w_vec )

          # print("prop_w_vec = ")
          # print(prop_w_vec)

          w_by_pred_sums_prop <-treepredmat %*% prop_w_vec  # rowSums(treepredmat %r*% prop_w_vec ) # rowSums(treepredmat %r+% prop_w_vec )

          y_hat <-  (1 - epsilon_w) * (y_hat_unnorm/att_weights_current_denoms) + epsilon_w*w_by_pred_sums
          y_hat_prop <-  (1 - epsilon_w) * (y_hat_unnorm/att_weights_current_denoms) + epsilon_w*w_by_pred_sums_prop

          # resids_orig <- y_scale - y_hat
          # resids_prop <- y_scale - y_hat_prop

          l_lik_orig <- sum(dnorm(x = y_scale,
                            mean = y_hat,
                            sd = sqrt(sigma2),
                            log = TRUE))

          l_lik_prop <- sum(dnorm(x = y_scale,
                            mean = y_hat_prop,
                            sd = sqrt(sigma2),
                            log = TRUE))



          w_mh_ratio <- exp( l_lik_prop  - l_lik_orig  +
                                 l_prior_prop - l_prior_orig +
                                 l_reverse_prob - l_prop_prob)


          if(is.na(w_mh_ratio)){

            print("v_w_vec = ")
            print(v_w_vec)

            print("step_size = ")
            print(step_size)

            print("i = ")
            print(i)

            print("j = ")
            print(j)

            print("v_j_prop = ")
            print(v_j_prop)

            print("prop_w_vec = ")
            print(prop_w_vec)

            print("l_lik_prop = ")
            print(l_lik_prop)

            print("l_lik_orig = ")
            print(l_lik_orig)

            print("l_prior_prop = ")
            print(l_prior_prop)

            print("l_prior_orig = ")
            print(l_prior_orig)

            print("l_reverse_prob = ")
            print(l_reverse_prob)

            print(" v_w_vec[j] = ")
            print( v_w_vec[j])

            print("reverse_step_size = ")
            print(reverse_step_size)

            print("l_prop_prob = ")
            print(l_prop_prob)


            print("l_prop_prob = ")
            print(l_prop_prob)


            print("a_w_vec[j] = ")
            print(a_w_vec[j])

            print("v_j_prop = ")
            print(v_j_prop)


            print("v_w_vec[j] = ")
            print(v_w_vec[j])


            # print("reverse_step_size = ")
            # print(reverse_step_size)
            #
            # print("l_prop_prob = ")
            # print(l_prop_prob)
            #
            #
            # print("l_prop_prob = ")
            # print(l_prop_prob)


          }


          if(runif(1) < w_mh_ratio){
            w_vec <- prop_w_vec
            v_w_vec <- prop_v_w_vec
          }

        } # end loop over j

        w_by_pred_sums <- treepredmat %*% w_vec  # rowSums(treepredmat %r*% w_vec ) # rowSums(treepredmat %r+% w_vec )


      } # end if dirichlet prior


    } # end w update




    # If at the right place, store everything
    if ((i > n_burn) & ((i - n_burn) %% n_thin) == 0) {
      curr <- (i - n_burn) / n_thin

      tree_store[[curr]] <- curr_trees
      sigma2_store[curr] <- sigma2
      y_hat_store[curr, ] <- y_hat
      var_count_store[curr, ] <- var_count
      s_prob_store[curr, ] <- s
      att_weights_store[[curr]] <- att_weights_current
      tau_store[curr] <- tau

      #store w_vec and tau_w
      if(include_w){
        wvec_store[curr,] <- w_vec
        tau_w_store[curr] <- tau_w
      }


    }
    # pb$tick()
  } # End loop through MCMC iterations

  cat("\n") # Make sure progress bar ends on a new line


  list_to_return <- list(
    trees = tree_store,
    sigma2 = sigma2_store * y_sd^2,
    y_hat =  (y_hat_store+(y_max + y_min)/2)*y_sd + y_mean,
    var_count_store = var_count_store,
    s = s_prob_store,
    # center = X_center,
    # scale = X_scale,
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
    # tau = tau,
    y_max = y_max,
    y_min = y_min,
    const_tree_weights = const_tree_weights,
    sq_num_features = sq_num_features,
    splitprob_as_weights = splitprob_as_weights,
    scale_x_funcs = scale_x_funcs,
    tau_store = tau_store,
    include_w = include_w
  )

  if(include_w){
    list_to_return$w_vecs <- wvec_store
    list_to_return$tau_w_store <- tau_w_store
    list_to_return$epsilon_w <- epsilon_w


  }

  return(list_to_return)
}
