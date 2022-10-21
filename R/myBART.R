#' @import Rcpp
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom rmutil 'ddoublepois'
#' @useDynLib AttBART, .registration = TRUE
#' @export

attbart = function(xtrain,
                   y,
                   sparse = TRUE,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1,
                   penalise_num_cov = TRUE,
                   lambda_cov = 0.4,
                   nu_cov = 2) {

  # NORMLIZE X
  # LOOK UP DBARTS, SoftBART, OR OTHER PACKAGES
  # MAYBE USE BART-IS NORMALIZATION

  x <- scale(xtrain)

  tempcenter <- attr(x, 'scaled:scale')
  tempscale <- attr(x, 'scaled:center')


  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {

      y_hat <- get_predictions(curr_trees, x, single_tree = FALSE)


      # print("y_hat = ")
      # print(y_hat)

      curr = (i - nburn)/nthin

      # print("curr = ")
      # print(curr)
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
    }

      # Start looping through trees
      for (j in 1:ntrees) {

        new_trees <- curr_trees

        # current_partial_residuals = y_scale - y_hat + tree_fits_store[,j]

        # Propose a new tree via grow/change/prune/swap
        # type = sample(c('grow', 'prune', 'change', 'swap'), 1)
        type = sample(c('grow', 'prune', 'change'), 1)
        if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

        # Generate a new tree based on the current
        new_trees[[j]] = update_tree(y = y_scale,
                                     X = x,
                                     type = type,
                                     curr_tree = curr_trees[[j]],
                                     node_min_size = node_min_size,
                                     s = s)

        # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior

        # Must account for attention weights given the current/old tree

        # to obtain current tree attention weight, need attention weights of all trees because normalized in softmax

        Aweights <- get_attention(curr_trees, x)

        # might need to obtain all tree predictions

        # matrix of tree predictions
        treepreds <- matrix(NA,
                            nrow = nrow(x),
                            ncol = length(curr_trees))

        for(tree_ind in 1:length(curr_trees)){

          # save unweighted predictions for tree
          treepreds[, tree_ind] <- get_predictions(curr_trees[[tree_ind]], x, single_tree = TRUE)

        }


        # create partial residuals conditional on the old tree

        #transformed current_partial_residuals
        # divide each partial residual by relevant attention weight of current/old tree
        # norm_partial_residuals <- current_partial_residuals/Aweights[,j]

        temp_par_resid <- y_scale - get_predictions(curr_trees, x, single_tree = FALSE) + treepreds[, j]*Aweights[,j]

        # if(!all(temp_par_resid == current_partial_residuals)){
        #   stop("Partial residuals calculation possibly incorrect")
        # }

        temp_par_resid <- temp_par_resid / Aweights[,j]


        # check that current_partial_residuals is correct


        #full conditional calculations must be adjusted to account for recaled errors

        # l_old = tree_full_conditional(curr_trees[[j]],
        #                               current_partial_residuals,
        #                               sigma2,
        #                               sigma2_mu) +
        #   get_tree_prior(curr_trees[[j]], alpha, beta) +
        #   get_num_cov_prior(curr_trees[[j]], lambda_cov, nu_cov, penalise_num_cov)

        # print("temp_par_resid = ")
        # print(temp_par_resid)

        l_old = att_tree_full_conditional(curr_trees[[j]],
                                          temp_par_resid,#current_partial_residuals,
                                          sigma2,
                                          sigma2_mu,
                                          Aweights[,j]) +
          get_tree_prior(curr_trees[[j]], alpha, beta) +
          get_num_cov_prior(curr_trees[[j]], lambda_cov, nu_cov, penalise_num_cov)

        # Must account for attention weights given the new tree



        # obtain attention weights for set of trees with nre tree

        Aweights_new <- get_attention(new_trees, x)

        # print("line 182 Aweights_new = ")
        # print(Aweights_new)
        #
        # print("line 185. get_predictions(new_trees, x, single_tree = FALSE) = ")
        # print(get_predictions(new_trees, x, single_tree = FALSE))
        #
        # print("line 185. get_predictions(new_trees[[j]], x, single_tree = TRUE)*Aweights_new[,j] = ")
        # print(get_predictions(new_trees[[j]], x, single_tree = TRUE)*Aweights_new[,j])

        # temp_par_resid_new <- y_scale - get_predictions(new_trees, x, single_tree = FALSE) +
        #   get_predictions(new_trees[[j]], x, single_tree = TRUE)*Aweights_new[,j]


        temp_par_resid_new <- y_scale - get_predictions_drop(new_trees, x, single_tree = FALSE,j)


        # print("line 184. temp_par_resid_new = ")
        # print(temp_par_resid_new)

        temp_par_resid_new <- temp_par_resid_new / Aweights_new[,j]

        # print("temp_par_resid_new = ")
        # print(temp_par_resid_new)


        # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_new = att_tree_full_conditional(new_trees[[j]],
                                          temp_par_resid_new,
                                          sigma2,
                                          sigma2_mu,
                                          Aweights_new[,j]) +
          get_tree_prior(new_trees[[j]], alpha, beta)  +
          get_num_cov_prior(new_trees[[j]], lambda_cov, nu_cov, penalise_num_cov)





        # Exponentiate the results above
        a = exp(l_new - l_old)


        if((i %%100) ==0){
          print("l_new = ")
          print(l_new)

          print("l_old = ")
          print(l_old)

          print("a = ")
          print(a)

        }



        temp_weights <- Aweights[,j]

        if(a > runif(1)) {
          curr_trees[[j]] = new_trees[[j]]
          temp_par_resid <- temp_par_resid_new
          temp_weights <- Aweights_new[,j]

          if (type =='change'){
            var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] - 1
            var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] + 1
          }

          if (type=='grow'){
            var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1 } # -1 because of the intercept in X

          if (type=='prune'){
            var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1 } # -1 because of the intercept in X
        }



        # NEED NEW DRAW OF TERMINAL NODE PARAMETERS



        # Update mu whether tree accepted or not
        curr_trees[[j]] = att_simulate_mu(curr_trees[[j]],
                                          temp_par_resid,#current_partial_residuals,
                                          sigma2,
                                          sigma2_mu,
                                          temp_weights)

      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      # y_hat = y_hat - tree_fits_store[,j] # subtract the old fit
      # y_hat = y_hat + current_fit # add the new fit

      #need to use all trees because weights depend on all trees
      y_hat <- get_predictions(curr_trees, x, single_tree = FALSE)

      tree_fits_store[,j] = current_fit # update the new fit

      } # End loop through trees

    # y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)
    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              var_count_store = var_count_store,
              s = s_prob_store,
              center = tempcenter,
              scale = tempscale,
              scaledtrainingdata = x
              ))

} # End main function



#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom rmutil 'ddoublepois'

cl_bart = function(x,
                   y,
                   sparse = TRUE,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1,
                   penalise_num_cov = TRUE,
                   lambda_cov = 0.4,
                   nu_cov = 2) {

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  n = length(y)
  p = ncol(x)
  s = rep(1/p, p)

  # Initial values
  z = ifelse(y == 0, -3, 3)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = z,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = pnorm(y_hat)
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      current_partial_residuals = z - y_hat + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      # type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      type = sample(c('grow', 'prune', 'change'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = z,
                                   X = x,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(curr_trees[[j]], alpha, beta) +
        get_num_cov_prior(curr_trees[[j]], lambda_cov, nu_cov, penalise_num_cov)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(new_trees[[j]], alpha, beta)  +
        get_num_cov_prior(new_trees[[j]], lambda_cov, nu_cov, penalise_num_cov)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] - 1
          var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] + 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1 } # -1 because of the intercept in X

        if (type=='prune'){
          var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1 } # -1 because of the intercept in X
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      y_hat = y_hat - tree_fits_store[,j] # subtract the old fit
      y_hat = y_hat + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Update z (latent variable)
    z = update_z(y, y_hat)

    # y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)
    sum_of_squares = sum((y - pnorm(y_hat))^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = n, nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store,
              y_hat = y_hat_store,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              var_count_store = var_count_store,
              s = s_prob_store
  ))

} # End main function
