# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to [...]     #
# -------------------------------------------------------------------------#
#
# 1. fill_tree_details: takes a tree matrix and returns the number of obs in
#       each node in it and the indices of each observation in each terminal node
# 2. get_predictions_no_w: gets the predicted values from a current set of trees
#       without the implementation of the attention parameters w
# 3. get_attention_no_w: calculates the attention weights from a current set of tres
#       without the implementation of the attention parameters w
# 4. get_children: takes a node and, if the node is terminal, returns the node.
#       If not, returns the children and calls the function again on the children
# 5. resample: an auxiliar function



# Fill_tree_details -------------------------------------------------------

fill_tree_details <- function(curr_tree, X) {
  # Collect right bits of tree
  tree_matrix <- curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix <- tree_matrix

  # Start with dummy node indices
  node_indices <- rep(1, nrow(X))

  # Exclude stump trees
  if (nrow(tree_matrix) > 1) {
    # For all but the top row (root node), find the number of observations falling into each one
    for (i in 2:nrow(tree_matrix)) {
      # Get the parent
      curr_parent <- tree_matrix[i, "parent"]

      # Find the split variable and value of the parent
      split_var <- tree_matrix[curr_parent, "split_variable"]
      split_val <- tree_matrix[curr_parent, "split_value"]

      # Find whether it's a left or right terminal node
      # left_or_right <- ifelse(tree_matrix[curr_parent, "child_left"] == i,
      #                         "left", "right"
      # )
      # if (left_or_right == "left") {
      if (tree_matrix[curr_parent, "child_left"] == i) {
          # If left use less than condition
        new_tree_matrix[i, "node_size"] <- sum(X[node_indices == curr_parent, split_var] < split_val)
        node_indices[node_indices == curr_parent][X[node_indices == curr_parent, split_var] < split_val] <- i
      } else {
        # If right use greater than condition
        # new_tree_matrix[i, "node_size"] <- sum(X[node_indices == curr_parent, split_var] >= split_val)
        # node_indices[node_indices == curr_parent][X[node_indices == curr_parent, split_var] >= split_val] <- i
        # remaining indices must be greater than or equal
        new_tree_matrix[i, "node_size"] <- sum(node_indices == curr_parent)
        node_indices[node_indices == curr_parent] <- i
      }
    }
  }

  return(list(
    tree_matrix = new_tree_matrix,
    node_indices = node_indices
  ))
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size = 1), ...]

# Rewrites the sample function into a function that allows "sampling" from a
#   vector of length 1, instead op applying the build-in feature of sample()
#   for single input values
sample.vec <- function(x, ...) x[sample(length(x), ...)]

update_s <- function(var_count, p, alpha_s) {
  s_ <- rdirichlet(1, as.vector((alpha_s / p ) + var_count))
  return(s_)
}

get_number_distinct_cov <- function(tree) {
  # Select the rows that correspond to internal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  # Get the covariates that are used to define the splitting rules
  num_distinct_cov <- length(unique(tree$tree_matrix[which_terminal, "split_variable"]))

  return(num_distinct_cov)
}


sample_move = function(curr_tree, i, nburn, trans_prob){

  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'),  1, prob = trans_prob)
  }
  return(type)
}


# get_children ------------------------------------------------------------

get_children <- function(tree_mat, parent) {
  # Create a holder for the children
  all_children <- NULL
  if (as.numeric(tree_mat[parent, "terminal"]) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not, get the current children
    curr_child_left <- as.numeric(tree_mat[parent, "child_left"])
    curr_child_right <- as.numeric(tree_mat[parent, "child_right"])
    # Return the children and also the children of the children recursively
    return(c(
      all_children,
      get_children(tree_mat, curr_child_left),
      get_children(tree_mat, curr_child_right)
    ))
  }
}

# Get predictions ---------------------------------------------------------

get_prediction_no_w <- function(tree, X) {
  # This is B_j in the paper (given tree j)

  # Stump trees just have one prediction value
  if (nrow(tree$tree_matrix) == 1) {
    predictions <- rep(tree$tree_matrix[1, "mu"], nrow(X))
  } else {
    # Loop through the node indices to get predictions
    predictions <- rep(NA, nrow(X))
    curr_X_node_indices <- tree$node_indices
    unique_node_indices <- unique(curr_X_node_indices)

    for (i in unique_node_indices) {
      predictions[curr_X_node_indices == i] <- tree$tree_matrix[i, "mu"] # tree$tree_matrix[[i, "mu"]]
    }
  }

  # Error message
  if (any(is.na(predictions))) {
    stop("Error in get_predictions. NA prediction values")
  }

  return(predictions)
}


# Get predictions test ---------------------------------------------------------

get_predictions_no_w_test = function(trees, X, single_tree = FALSE, tau, feature_weighting, const_tree_weights,
                                     sq_num_features,
                                     splitprob_as_weights,
                                     s,
                                     test_binary = TRUE,
                                     include_w,
                                     epsilon_w,
                                     w_vec) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  if(const_tree_weights){
    att_weights <- matrix(1/length(trees), nrow = nrow(X), ncol = length(trees))
  }else{
    att_weights <- get_attention_no_w(trees, X, tau, feature_weighting, sq_num_features, splitprob_as_weights,
                                      s, test_binary)
  }

  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    predictions <- rep(0, nrow(X))
    if(include_w){
      w_by_preds <- rep(0, nrow(X))
    }
    for(m in 1:length(trees)){

      temptree <- trees[[m]]

      # Loop through the node indices to get predictions
      treepredictions = rep(NA, nrow(X))
      unique_node_indices = unique(temptree$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(temptree, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        treepredictions[curr_X_node_indices == unique_node_indices[i]] =
          temptree$tree_matrix[unique_node_indices[i], 'mu']
      }

      predictions <- predictions + treepredictions*att_weights[,m]
      if(include_w){
        w_by_preds <- w_by_preds + w_vec[m]*treepredictions
      }
    }

    if(include_w){
      predictions <- (1- epsilon_w)*predictions + epsilon_w*w_by_preds
    }

    #   # Do a recursive call to the function
    # partial_trees = trees
    # att_weight_vec <- att_weights_current[,1]
    # att_weights_current <- att_weights_current[,-1, drop = FALSE]
    # partial_trees[[1]] = NULL # Blank out that element of the list
    # predictions = get_predictions_no_w_test(trees[[1]], X, single_tree = TRUE)*att_weight_vec  +
    #   get_predictions_no_w_test(partial_trees, X,
    #                   single_tree = length(partial_trees) == 1)
    # #single_tree = !is.null(names(partial_trees)))
    # # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}



# get unnormalized attention for one tree --------------------------------------------------


get_unnorm_att_1tree_no_w <- function(curr_tree, X, tau = 1, feature_weighting, sq_num_features,
                               splitprob_as_weights,
                               s,
                               test_binary = TRUE) {
  # To calculate attention weights:
  #  1. For each tree, find the leaf for each training observation x_i
  #  2. For each tree, obtain leaf covariate means A_j(x_i) (vectors)
  #  3. For each tree and observation, (a) subtract the covariate mean of the relevant leaf,
  #        (b) square the difference, (c) divide by 2*tau and (d) then take the exponential
  #  4. After all the calculations for all trees and observations, normalise the weights of
  #        the observations over the trees to sum to 1 using the softmax function



  # Matrix of observation and tree-specific unnormalised weights
  # weight_exp_matrix <- matrix(NA,
  #                             nrow = nrow(X),
  #                             ncol = length(trees)
  # )

  # Get attention weights for each tree j
  # for (j in 1:length(trees)) {
  #
  #
  #   # Find the current terminal nodes (leaves)
  #   curr_tree <- trees[[j]]

  if(test_binary){
    curr_tree <- fill_tree_details(curr_tree, X)
  }else{
    # curr_tree <- trees[[j]]
  }

    if(!splitprob_as_weights){
      # Feature weighting
      if (feature_weighting) {
        tree_splitvars <- curr_tree$tree_matrix[, "split_variable"]
        if (any(!is.na(tree_splitvars))) {
          # Remove the NAs
          splitvars <- unique(na.omit(tree_splitvars))
          # splitvars <- sort(unique(na.omit(tree_splitvars)), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
        } else {
          # Else, return all variables
          splitvars <- 1:ncol(X)
        }
      } else {
        splitvars <- NA
      }
    }

    X_node_indices <- curr_tree$node_indices
    # unique_leaf_indices <- unique(X_node_indices)



    # Calculate the L2 norm and unnormalised weights splitvars is NA when no
    # feature weighting is specified and a vector of either length 1 up to
    # length p

    if(splitprob_as_weights){

      # # Initialise mean matrix
      # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
      #
      # # Calculate the means A_j(x_i)
      # for (l in unique_leaf_indices) {
      #   leaf_indices <- X_node_indices == l
      #   mean_matrix[leaf_indices, ] <- rep(colMeans(X[leaf_indices, ]), each = sum(leaf_indices))
      # }

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[X_node_indices,]

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[as.character(X_node_indices),]

      X_diff_sq <- (fmean(X, X_node_indices, TRA = "-") %r*% s )^2

      # X_diff_sq <- ( ((X - mean_matrix) %r*% s )   ^2)
      X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # if(sq_num_features){
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # }else{
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)
      # }

      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      # weight_exp_matrix[, j] <- -1 * X_L2_sq / (2 * tau)
      weight_exp_vec <- exp(-1 * X_L2_sq / (2 * tau) )

    }else{

      if (any(!is.na(splitvars))) {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = length(splitvars))
        # t_mean_matrix <- matrix(NA, nrow = length(splitvars) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, splitvars, drop = FALSE])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # colmeanmat <- fmean(X[, splitvars, drop = FALSE], X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]
        #
        # # X_diff_sq <- (X[, splitvars] - mean_matrix[, splitvars])^2
        # X_diff_sq <- (X[, splitvars, drop = FALSE] - mean_matrix)^2


        X_diff_sq <- fmean(X[, splitvars, drop = FALSE], X_node_indices, TRA = "-")^2


        # X_diff_sq <- (X[, splitvars, drop = FALSE] - t(t_mean_matrix))^2
        # If splitvars only has one variable, then rowSums does not work
        if (length(splitvars) == 1) {
          X_L2_sq <- X_diff_sq
        } else {
          if(sq_num_features){
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
          }else{
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
          }
        }
      } else {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
        # t_mean_matrix <- matrix(NA, nrow = ncol(X) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, ])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # mean_matrix <- fmean(X, X_node_indices)
        # colmeanmat <- fmean(X, X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]

        X_diff_sq <- fmean(X, X_node_indices, TRA = "-")^2

        # X_diff_sq <- ((X - mean_matrix)^2)
        # X_diff_sq <- ((X - t(t_mean_matrix))^2)
        if(sq_num_features){
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
        }else{
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
        }
      }

      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      weight_exp_vec <- exp(-1 * X_L2_sq / (2 * tau))

    }
  # }

  # weight_exp_matrix <- exp(weight_exp_matrix)

  # Normalize distances to obtain softmax weights by divide each row by the row sum
  # att_weights <- t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))
  # att_weights <- weight_exp_matrix/rowSums(weight_exp_matrix)  #t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))

  return(weight_exp_vec)
}




# Get unnormalized attention for all treesv-----------------------------------------------------------

get_unnorm_att_all_no_w <- function(trees, X, tau = 1, feature_weighting, sq_num_features,
                               splitprob_as_weights,
                               s,
                               test_binary = TRUE) {
  # To calculate attention weights:
  #  1. For each tree, find the leaf for each training observation x_i
  #  2. For each tree, obtain leaf covariate means A_j(x_i) (vectors)
  #  3. For each tree and observation, (a) subtract the covariate mean of the relevant leaf,
  #        (b) square the difference, (c) divide by 2*tau and (d) then take the exponential
  #  4. After all the calculations for all trees and observations, normalise the weights of
  #        the observations over the trees to sum to 1 using the softmax function



  # Matrix of observation and tree-specific unnormalised weights
  weight_exp_matrix <- matrix(NA,
                              nrow = nrow(X),
                              ncol = length(trees)
  )

  # Get attention weights for each tree j
  for (j in 1:length(trees)) {


    # Find the current terminal nodes (leaves)
    if(test_binary){
      curr_tree <- fill_tree_details(trees[[j]], X)
    }else{
      curr_tree <- trees[[j]]
    }

    if(!splitprob_as_weights){
      # Feature weighting
      if (feature_weighting) {
        tree_splitvars <- curr_tree$tree_matrix[, "split_variable"]
        if (any(!is.na(tree_splitvars))) {
          # Remove the NAs
          splitvars <- unique(na.omit(tree_splitvars))
          # splitvars <- sort(unique(na.omit(tree_splitvars)), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
        } else {
          # Else, return all variables
          splitvars <- 1:ncol(X)
        }
      } else {
        splitvars <- NA
      }
    }

    X_node_indices <- curr_tree$node_indices
    # unique_leaf_indices <- unique(X_node_indices)



    # Calculate the L2 norm and unnormalised weights splitvars is NA when no
    # feature weighting is specified and a vector of either length 1 up to
    # length p

    if(splitprob_as_weights){

      # # Initialise mean matrix
      # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
      #
      # # Calculate the means A_j(x_i)
      # for (l in unique_leaf_indices) {
      #   leaf_indices <- X_node_indices == l
      #   mean_matrix[leaf_indices, ] <- rep(colMeans(X[leaf_indices, ]), each = sum(leaf_indices))
      # }

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[X_node_indices,]

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[as.character(X_node_indices),]

      X_diff_sq <- (fmean(X, X_node_indices, TRA = "-") %r*% s )^2

      # X_diff_sq <- ( ((X - mean_matrix) %r*% s )   ^2)
      X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # if(sq_num_features){
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # }else{
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)
      # }

      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      weight_exp_matrix[, j] <- -1 * X_L2_sq / (2 * tau)

    }else{

      if (any(!is.na(splitvars))) {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = length(splitvars))
        # t_mean_matrix <- matrix(NA, nrow = length(splitvars) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, splitvars, drop = FALSE])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # colmeanmat <- fmean(X[, splitvars, drop = FALSE], X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]
        #
        # # X_diff_sq <- (X[, splitvars] - mean_matrix[, splitvars])^2
        # X_diff_sq <- (X[, splitvars, drop = FALSE] - mean_matrix)^2


        X_diff_sq <- fmean(X[, splitvars, drop = FALSE], X_node_indices, TRA = "-")^2


        # X_diff_sq <- (X[, splitvars, drop = FALSE] - t(t_mean_matrix))^2
        # If splitvars only has one variable, then rowSums does not work
        if (length(splitvars) == 1) {
          X_L2_sq <- X_diff_sq
        } else {
          if(sq_num_features){
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
          }else{
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
          }
        }
      } else {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
        # t_mean_matrix <- matrix(NA, nrow = ncol(X) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, ])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # mean_matrix <- fmean(X, X_node_indices)
        # colmeanmat <- fmean(X, X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]

        X_diff_sq <- fmean(X, X_node_indices, TRA = "-")^2

        # X_diff_sq <- ((X - mean_matrix)^2)
        # X_diff_sq <- ((X - t(t_mean_matrix))^2)
        if(sq_num_features){
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
        }else{
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
        }
      }

      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      weight_exp_matrix[, j] <- -1 * X_L2_sq / (2 * tau)

    }
  }

  weight_exp_matrix <- exp(weight_exp_matrix)

  # Normalize distances to obtain softmax weights by divide each row by the row sum
  # att_weights <- t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))
  att_weights <- weight_exp_matrix#/rowSums(weight_exp_matrix)  #t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))

  return(att_weights)
}

# Get attention -----------------------------------------------------------

get_attention_no_w <- function(trees, X, tau = 1, feature_weighting, sq_num_features,
                               splitprob_as_weights,
                               s,
                               test_binary = TRUE) {
  # To calculate attention weights:
  #  1. For each tree, find the leaf for each training observation x_i
  #  2. For each tree, obtain leaf covariate means A_j(x_i) (vectors)
  #  3. For each tree and observation, (a) subtract the covariate mean of the relevant leaf,
  #        (b) square the difference, (c) divide by 2*tau and (d) then take the exponential
  #  4. After all the calculations for all trees and observations, normalise the weights of
  #        the observations over the trees to sum to 1 using the softmax function



  # Matrix of observation and tree-specific unnormalised weights
  weight_exp_matrix <- matrix(NA,
                              nrow = nrow(X),
                              ncol = length(trees)
  )

  # Get attention weights for each tree j
  for (j in 1:length(trees)) {


    # Find the current terminal nodes (leaves)
    if(test_binary){
      curr_tree <- fill_tree_details(trees[[j]], X)
    }else{
      curr_tree <- trees[[j]]
    }

    if(!splitprob_as_weights){
      # Feature weighting
      if (feature_weighting) {
        tree_splitvars <- curr_tree$tree_matrix[, "split_variable"]
        if (any(!is.na(tree_splitvars))) {
          # Remove the NAs
          splitvars <- unique(na.omit(tree_splitvars))
          # splitvars <- sort(unique(na.omit(tree_splitvars)), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
          # splitvars <- sort(na.omit(tree_splitvars), na.last = TRUE)
        } else {
          # Else, return all variables
          splitvars <- 1:ncol(X)
        }
      } else {
        splitvars <- NA
      }
    }

    X_node_indices <- curr_tree$node_indices
    # unique_leaf_indices <- unique(X_node_indices)



    # Calculate the L2 norm and unnormalised weights splitvars is NA when no
    # feature weighting is specified and a vector of either length 1 up to
    # length p

    if(splitprob_as_weights){

      # # Initialise mean matrix
      # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
      #
      # # Calculate the means A_j(x_i)
      # for (l in unique_leaf_indices) {
      #   leaf_indices <- X_node_indices == l
      #   mean_matrix[leaf_indices, ] <- rep(colMeans(X[leaf_indices, ]), each = sum(leaf_indices))
      # }

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[X_node_indices,]

      # colmeanmat <- fmean(X, X_node_indices)
      # mean_matrix <- colmeanmat[as.character(X_node_indices),]

      X_diff_sq <- (fmean(X, X_node_indices, TRA = "-") %r*% s )^2

      # X_diff_sq <- ( ((X - mean_matrix) %r*% s )   ^2)
      X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # if(sq_num_features){
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)^2
      # }else{
      #   X_L2_sq <- rowSums(X_diff_sq)#/ length(splitvars)
      # }

      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      weight_exp_matrix[, j] <- -1 * X_L2_sq / (2 * tau)

    }else{

      if (any(!is.na(splitvars))) {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = length(splitvars))
        # t_mean_matrix <- matrix(NA, nrow = length(splitvars) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, splitvars, drop = FALSE])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # colmeanmat <- fmean(X[, splitvars, drop = FALSE], X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]
        #
        # # X_diff_sq <- (X[, splitvars] - mean_matrix[, splitvars])^2
        # X_diff_sq <- (X[, splitvars, drop = FALSE] - mean_matrix)^2


        X_diff_sq <- fmean(X[, splitvars, drop = FALSE], X_node_indices, TRA = "-")^2


        # X_diff_sq <- (X[, splitvars, drop = FALSE] - t(t_mean_matrix))^2
        # If splitvars only has one variable, then rowSums does not work
        if (length(splitvars) == 1) {
          X_L2_sq <- X_diff_sq
        } else {
          if(sq_num_features){
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
          }else{
            X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
          }
        }
      } else {

        # # Initialise mean matrix
        # # mean_matrix <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
        # t_mean_matrix <- matrix(NA, nrow = ncol(X) , ncol = nrow(X))
        #
        # # Calculate the means A_j(x_i)
        # for (l in unique_leaf_indices) {
        #   leaf_indices <- X_node_indices == l
        #   tempcolmeans <- colMeans(X[leaf_indices, ])
        #   t_mean_matrix[, leaf_indices] <- tempcolmeans
        #
        #   # mean_matrix[leaf_indices, ] <- rep(tempcolmeans, each = sum(leaf_indices))
        # }

        # mean_matrix <- fmean(X, X_node_indices)
        # colmeanmat <- fmean(X, X_node_indices)
        # mean_matrix <- colmeanmat[as.character(X_node_indices),]

        X_diff_sq <- fmean(X, X_node_indices, TRA = "-")^2

        # X_diff_sq <- ((X - mean_matrix)^2)
        # X_diff_sq <- ((X - t(t_mean_matrix))^2)
        if(sq_num_features){
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)^2
        }else{
          X_L2_sq <- rowSums(X_diff_sq)/ length(splitvars)
        }
      }


      # weight_exp_matrix[, j] <- exp(-1 * X_L2_sq / (2 * tau))
      weight_exp_matrix[, j] <- -1 * X_L2_sq / (2 * tau)

    }
  }

  weight_exp_matrix <- exp(weight_exp_matrix)

  # Normalize distances to obtain softmax weights by divide each row by the row sum
  # att_weights <- t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))
  att_weights <- weight_exp_matrix/rowSums(weight_exp_matrix)  #t(apply(weight_exp_matrix, 1, function(x) x / sum(x)))

  return(att_weights)
}

# MH acceptance probability -----------------------------------------------

# This function returns the number of second generation nodes "w"
get_gen2 <- function(tree) {
  if (nrow(tree$tree_matrix) == 1) {
    w <- 0
  } else {
    # indeces <- which(tree$tree_matrix[, "terminal"] == 1)
    # Determine the parent for each terminal node and find the duplicated parents
    # w <- as.numeric(sum(duplicated(tree$tree_matrix[indeces,'parent'])))
    # parents <- tree$tree_matrix[indeces, "parent"]
    parents <- tree$tree_matrix[tree$tree_matrix[, "terminal"] == 1, "parent"]
    w <- parents[duplicated(parents)]
  }
  return(w)
}

# This function returns a vector of the number of unique observations for each
# covariate in the given terminal node
get_n_unique <- function(tree, X, node) {
  p <- ncol(X)
  n_unique <- sapply(1:p, function(q) {
    length(unique(X[tree$node_indices == node, q]))
  })
  return(n_unique)
}

# This function returns the depth of a given node in the tree, given depth=0 for
# the root node
get_depth <- function(tree, node) {
  depth_chosen_node <- 0
  current_row <- node
  while (current_row > 1) {
    current_row <- tree$tree_matrix[current_row, "parent"]
    depth_chosen_node <- depth_chosen_node + 1
  }
  return(depth_chosen_node)
}

# This function finds and returns the node that was pruned in the current tree
get_pruned_node <- function(curr_tree, pruned_tree) {
  # The chosen node that it used to prune is the first row where the new and current tree
  #   matrices do not agree, easiest to search through the column "split_variable"
  pruned_node <- 1
  i <- 1
  while (is.na(curr_tree$tree_matrix[i, "split_variable"]) == is.na(pruned_tree$tree_matrix[i, "split_variable"])) {
    i <- i + 1
    pruned_node <- pruned_node + 1
  }
  return(pruned_node)
}

# This function returns the Metropolis-Hastings acceptance probability in line
# with Kapelner, A., and Bleich, J. (2013). "bartMachine: Machine learning with
# Bayesian additive regression trees."

get_MH_probability <- function(X, curr_tree, new_tree,
                               att_weights_current, att_weights_new,
                               curr_partial_resid_rescaled,
                               new_partial_resid_rescaled,
                               type, trans_prob,
                               alpha, beta,
                               mu_mu, sigma2_mu, sigma2,
                               node_min_size) {

  # Number of terminal nodes in current tree
  b_j <- sum(curr_tree$tree_matrix[, "terminal"])

  # Get the tree type probabilities
  # if(nrow(curr_tree$tree_matrix) == 1 ){
  #   prob_grow <-  1 #trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }else{
  #   prob_grow <- trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }

  prob_grow <- trans_prob[1]
  prob_prune <- trans_prob[2]

  # print("type = ")
  # print(type)
  # Type is either "grow", "prune" or "change"
  if (type == "grow") {
    # The chosen node that it used to grow is the parent of the newly added last two rows in the NEW tree
    node_to_split <- as.numeric(new_tree$tree_matrix[nrow(new_tree$tree_matrix), "parent"])
    # Get the number of unique values for each covariate given the chosen node to split on in the CURRENT tree
    # n_unique <- get_n_unique(curr_tree, X, node_to_split)
    # Get the chosen covariate that it used to split on, which is the split variable in the NEW tree
    # covariate_to_split <- new_tree$tree_matrix[node_to_split, "split_variable"]
    # Obtain the additional parameters needed for the tree and transition ratios
    # p_j <- sum(n_unique >= 2 * node_min_size) # Number of available covariates to split on
    # n_j <- n_unique[covariate_to_split] - (node_min_size > 0) * (2 * node_min_size - 1) # Number of available values to split on given the chosen covariate, adjusted to the minimum leaf size
    w_2_star <- length(get_gen2(new_tree)) # Number of second generation nodes in the NEW tree

    # print("node_to_split = ")
    # print(node_to_split)
    #
    # print("curr_tree = ")
    # print(curr_tree)
    #
    # print("new_tree = ")
    # print(new_tree)

    depth_jl <- get_depth(curr_tree, node_to_split) # The depth of the chosen node to split on (in either the current or new tree)

    # tree_ratio <- alpha * (1 - alpha / ((2 + depth_jl)^beta))^2 / (((1 + depth_jl)^beta - alpha) * p_j * n_j)
    # transition_ratio <- prob_prune / prob_grow * b_j * p_j * n_j / w_2_star
    # print("depth_jl = ")
    # print(depth_jl)
    #
    # print("alpha * (1 - alpha / ((2 + depth_jl)^beta))^2 = ")
    # print(alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    #
    # print("((1 + depth_jl)^beta - alpha)  = ")
    # print(((1 + depth_jl)^beta - alpha) )
    #
    # print("alpha = ")
    # print(alpha)
    #
    # print("beta = ")
    # print(beta)

    tree_ratio <- alpha * (1 - (alpha / ((2 + depth_jl)^beta)) )^2 / (((1 + depth_jl)^beta - alpha) )
    transition_ratio <- (prob_prune / prob_grow) * (b_j  / w_2_star) # maybe this should be (w_2_star + 1) ?
  } else if (type == "prune") {
    # Finding the node used to prune is a bit more difficult then finding the node to split on
    pruned_node <- get_pruned_node(curr_tree, new_tree)
    # Get the number of unique values for each covariate given the chosen node it pruned on in the NEW tree
    # n_unique <- get_n_unique(new_tree, X, pruned_node)
    # Get the covariate that was used in the pruned node, which is the split variable in the CURRENT tree
    # pruned_covariate <- curr_tree$tree_matrix[pruned_node, "split_variable"]
    # Obtain the additional parameters needed for the tree and transition ratios
    # p_j_star <- sum(n_unique >= 2 * node_min_size) # Number of available covariates to split on
    # n_j_star <- n_unique[pruned_covariate] - (node_min_size > 0) * (2 * node_min_size - 1) # Number of available values to split on given the chosen covariate, adjusted to the minimum leaf size
    w_2 <- length(get_gen2(curr_tree)) # Number of second generation nodes in the CURRENT tree

    # print("pruned_node = ")
    # print(pruned_node)
    #
    # print("new_tree = ")
    # print(new_tree)

    depth_jl <- get_depth(new_tree, pruned_node) # The depth of the chosen node to split on (in either the current or new tree)

    # tree_ratio <- ((1 + depth_jl)^beta - alpha) * p_j_star * n_j_star / (alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    # transition_ratio <- prob_grow / prob_prune * w_2 / ((b_j - 1) * p_j_star * n_j_star)

    tree_ratio <- ((1 + depth_jl)^beta - alpha)  / (alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    transition_ratio <- (prob_grow / prob_prune) * (w_2 / ((b_j - 1) ))
  } else {
    # For the changing step, the tree and transition ratios cancel out into a factor of 1
    transition_ratio <- 1
    tree_ratio <- 1
  }
  l_new <- get_logL(new_tree, new_partial_resid_rescaled, att_weights_new, mu_mu, sigma2_mu, sigma2)
  l_old <- get_logL(curr_tree, curr_partial_resid_rescaled, att_weights_current, mu_mu, sigma2_mu, sigma2)
  r <- exp(l_new - l_old) * transition_ratio * tree_ratio

  # print("l_new = ")
  # print(l_new)
  #
  # print("l_old = ")
  # print(l_old)
  #
  # print("transition_ratio = ")
  # print(transition_ratio)
  #
  # print("tree_ratio = ")
  # print(tree_ratio)

  return(min(1, r))
}



# This function returns the Metropolis-Hastings acceptance probability in line
# with Kapelner, A., and Bleich, J. (2013). "bartMachine: Machine learning with
# Bayesian additive regression trees."
get_MH_probability2 <- function(X, curr_tree, new_tree,
                               att_weights_current, att_weights_new,
                               curr_partial_resid_rescaled,
                               new_partial_resid_rescaled,
                               type, trans_prob,
                               alpha, beta,
                               mu_mu, sigma2_mu, sigma2,
                               node_min_size) {
  # Number of terminal nodes in current tree
  # b_j <- sum(curr_tree$tree_matrix[, "terminal"])

  # Get the tree type probabilities
  # if(nrow(curr_tree$tree_matrix) == 1 ){
  #   prob_grow <-  1 #trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }else{
  #   prob_grow <- trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }

  prob_grow <- trans_prob[1]
  prob_prune <- trans_prob[2]
  l_new <- get_logL(new_tree, new_partial_resid_rescaled, att_weights_new, mu_mu, sigma2_mu, sigma2)
  l_old <- get_logL(curr_tree, curr_partial_resid_rescaled, att_weights_current, mu_mu, sigma2_mu, sigma2)

  if(type == 'grow'){
    # a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) )*
    #   ratio_grow(curr_tree, new_tree) * (prob_prune / prob_grow)
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log(ratio_grow(curr_tree, new_tree ) ))*
       (prob_prune / prob_grow)
  } else if(type == 'prune'){
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log( ratio_prune(curr_tree, new_tree) ))*
     (prob_grow / prob_prune)
  } else{
    a = exp(l_new - l_old)
  }

  if(is.na(a)){
    print("get_tree_prior(new_tree, alpha, beta) = ")
    print(get_tree_prior(new_tree, alpha, beta))

    print("get_tree_prior(curr_tree, alpha, beta) ) = ")
    print(get_tree_prior(curr_tree, alpha, beta) )


    print("ratio_grow(curr_tree, new_tree ) ) = ")
    print(ratio_grow(curr_tree, new_tree ) )

    print("ratio_prune(curr_tree, new_tree)  = ")
    print(ratio_prune(curr_tree, new_tree) )

    print("new_tree = ")
    print(new_tree)

    print("curr_tree = ")
    print(curr_tree)

    print("alpha = ")
    print(alpha)

    print("beta = ")
    print(beta)

    print("l_new = ")
    print(l_new)

    print("l_old = ")
    print(l_old)


  }

  # r <- exp(l_new - l_old) * transition_ratio * tree_ratio
  return(min(1, a))
}

get_logL <- function(tree, residuals, att_weights, mu_mu, sigma2_mu, sigma2) {
  terminal_nodes <- which(tree$tree_matrix[, "terminal"] == 1)
  # logL <- 0

  nj = tree$tree_matrix[terminal_nodes,'node_size']


  att_weights_sq = att_weights^2
  sum_1_vec = fsum(att_weights_sq, tree$node_indices)
  sum_2_vec = fsum(att_weights_sq * residuals, tree$node_indices)
  # sum_2_vec = fsum(att_weights_sq , tree$node_indices, residuals)
  logsum_attw_vec = fsum(log(att_weights), tree$node_indices)
  awsq_by_resid2_vec = fsum(att_weights_sq * residuals^2, tree$node_indices)
  # awsq_by_resid2_vec = fsum(att_weights_sq , tree$node_indices, residuals^2)

  # print("terminal_nodes = ")
  # print(terminal_nodes)
  # # for (l in terminal_nodes) {
  # for (l in 1:length(terminal_nodes)) {
  #   # node_indices_l <- which(tree$node_indices == l)
  #   # att_weights_l <- att_weights[node_indices_l]
  #   # residuals_l <- residuals[node_indices_l]
  #   #
  #   # # Pre-calculations
  #   # att_weights_l_sq <- att_weights_l^2
  #   # sum_1 <- sum(att_weights_l_sq)
  #   # sum_2 <- sum(att_weights_l_sq * residuals_l)
  #
  #   # Calculate the log-likelihood of node l and add it to the total. The terms
  #   # of 0.5 * log(sigma2) and mu_mu^2 / (2 * sigma2_mu) can be removed as these
  #   # cancel out, see the full logL in the commented lines below
  #   # logL <- logL - 0.5 * length(node_indices_l) * log(2 * sigma2 * pi) + sum(log(att_weights_l)) -
  #   #   0.5 * log(sigma2_mu * sum_1 + sigma2) - sum(att_weights_l_sq * residuals_l^2) / (2 * sigma2) +
  #   #   ((sigma2_mu / sigma2) * sum_2^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2) / (2 * (sigma2_mu * sum_1 + sigma2))
  #
  #   # logL <- logL - 0.5 * length(node_indices_l) * log(2 * sigma2 * pi) + sum(log(att_weights_l)) +
  #   #   0.5 * log(sigma2) - 0.5 * log(sigma2_mu * sum(att_weights_l_sq) + sigma2) - sum(att_weights_l_sq * residuals_l^2) / (2 * sigma2) -
  #   #   mu_mu^2 / (2 * sigma2_mu) + ((sigma2_mu / sigma2) * sum(att_weights_l_sq * residuals_l)^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum(att_weights_l_sq * residuals_l)) / (2 * (sigma2_mu * sum(att_weights_l_sq) + sigma2))
  #   # tempval = - 0.5 * nj[l] * log(2 * sigma2 * pi) + logsum_attw_vec[l] -
  #   #   0.5 * log(sigma2_mu * sum_1_vec[l] + sigma2) - awsq_by_resid2_vec[l] / (2 * sigma2) +
  #   #   ((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]) / (2 * (sigma2_mu * sum_1_vec[l] + sigma2))
  #
  #   # print("tempval = ")
  #   # print(tempval)
  #   #
  #   # print("logL = ")
  #   # print(logL)
  #   logL <- logL - 0.5 * nj[l] * log(2 * sigma2 * pi) + logsum_attw_vec[l] -
  #     0.5 * log(sigma2_mu * sum_1_vec[l] + sigma2) - awsq_by_resid2_vec[l] / (2 * sigma2) +
  #     ((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]) / (2 * (sigma2_mu * sum_1_vec[l] + sigma2))
  #
  #   if(is.na(logL)){
  #     # print("nj[l] = ")
  #     # print(nj[l])
  #     # print("logsum_attw_vec[l] = ")
  #     # print(logsum_attw_vec[l])
  #     # print("sum_1_vec[l] = ")
  #     # print(sum_1_vec[l])
  #     # print("sum_2_vec[l] = ")
  #     # print(sum_2_vec[l])
  #     # print("awsq_by_resid2_vec[l] = ")
  #     # print(awsq_by_resid2_vec[l])
  #     #
  #     #
  #     # print("0.5 * nj[l] * log(2 * sigma2 * pi) = ")
  #     # print(0.5 * nj[l] * log(2 * sigma2 * pi))
  #     #
  #     # print("0.5 * log(sigma2_mu * sum_1_vec[l] + sigma2)  = ")
  #     # print(0.5 * log(sigma2_mu * sum_1_vec[l] + sigma2) )
  #     #
  #     # print("awsq_by_resid2_vec[l] / (2 * sigma2) = ")
  #     # print(awsq_by_resid2_vec[l] / (2 * sigma2))
  #     #
  #     # print("((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]) / (2 * (sigma2_mu * sum_1_vec[l] + sigma2)) = ")
  #     # print(((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]) / (2 * (sigma2_mu * sum_1_vec[l] + sigma2)))
  #     #
  #     # print("((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]) = ")
  #     # print(((sigma2_mu / sigma2) * sum_2_vec[l]^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec[l]))
  #     #
  #     # print("2 * (sigma2_mu * sum_1_vec[l] + sigma2) = ")
  #     # print(2 * (sigma2_mu * sum_1_vec[l] + sigma2))
  #     #
  #     # print("sum_2_vec[l]^2 = ")
  #     # print(sum_2_vec[l]^2)
  #     #
  #     # # print("tempval = ")
  #     # # print(tempval)
  #     #
  #     # print("l = ")
  #     # print(l)
  #
  #     stop("NA logL")
  #     # print("logL = ")
  #     # print(logL)
  #     # print("logL = ")
  #     # print(logL)
  #   }
  #
  # }

  logL <-  sum( - 0.5 * nj * log(2 * sigma2 * pi) + logsum_attw_vec -
    0.5 * log(sigma2_mu * sum_1_vec + sigma2) - awsq_by_resid2_vec / (2 * sigma2) +
    ((sigma2_mu / sigma2) * sum_2_vec^2 + (sigma2 / sigma2_mu) * mu_mu^2 + 2 * mu_mu * sum_2_vec) / (2 * (sigma2_mu * sum_1_vec + sigma2))
    )

  # print("logL = ")
  # print(logL)
  # print("logL2 = ")
  # print(logL2)

  # if(logL == logL2){
  #   print("logLs equal :)")
  # }
  # if(logL != logL2){
  #   print("logLs not equal ")
  #
  #   print("logL = ")
  #   print(logL)
  #   print("logL2 = ")
  #   print(logL2)
  # }

  return(logL)
}


# Update node parameters --------------------------------------------------

att_simulate_mu <- function(tree, residuals, att_weights, mu_mu, sigma2_mu, sigma2) {
  tree$tree_matrix[, "mu"] <- NA

  terminal_nodes <- which(tree$tree_matrix[, "terminal"] == 1)

  nj = tree$tree_matrix[terminal_nodes,'node_size']


  att_weights_sq = att_weights^2
  sum_1_vec = fsum(att_weights_sq, tree$node_indices)
  sum_2_vec = fsum(att_weights_sq * residuals, tree$node_indices)
  # sum_2_vec = fsum(att_weights_sq , tree$node_indices, residuals)
  # logsum_attw_vec = fsum(log(att_weights), tree$node_indices)
  # awsq_by_resid2_vec = fsum(att_weights_sq * residuals^2, tree$node_indices)
  # awsq_by_resid2_vec = fsum(att_weights_sq , tree$node_indices, residuals^2)


  # for (l in terminal_nodes) {
  #   # Get leaf weights and residuals
  #   node_indices_l <- which(tree$node_indices == l)
  #   att_weights_l_sq <- att_weights[node_indices_l]^2
  #   residuals_l <- residuals[node_indices_l]
  #
  #   # Sample mu
  #   var <- sigma2 * sigma2_mu / (sigma2 + sigma2_mu * sum(att_weights_l_sq))
  #   mean <- var * (sigma2 * mu_mu + sigma2_mu * sum(att_weights_l_sq * residuals_l)) / (sigma2_mu * sigma2)
  #
  #   mu <- rnorm(1, mean = mean, sd = sqrt(var))
  #
  #   # Update mu
  #   tree$tree_matrix[l, "mu"] <- mu
  # }


  # tree$tree_matrix[,'mu'] <- NA

  var_vec <- sigma2 * sigma2_mu / (sigma2 + sigma2_mu * sum_1_vec)
  mean_vec <- var_vec * (sigma2 * mu_mu + sigma2_mu * sum_2_vec) / (sigma2_mu * sigma2)

  tree$tree_matrix[terminal_nodes, "mu"] <- rnorm(length(terminal_nodes),
                                                  mean = mean_vec,
                                                  sd = sqrt(var_vec))

  # if (any(is.na(tree$tree_matrix[terminal_nodes, "mu"]))) {
  #   stop("Some terminal nodes do not have a prediction!")
  # }

  return(tree)
}

update_sigma2 <- function(S, n, nu, lambda) {
  sigma2 <- 1 / rgamma(1, shape = (n + nu) / 2, rate = (S + nu * lambda) / 2)
  return(sigma2)
}




update_alpha <- function(s, alpha_scale, alpha_a, alpha_b) {

  # create inputs for likelihood

  log_s <- log(s)
  mean_log_s <- mean(log_s)
  p <- length(s)
  # alpha_scale   # denoted by lambda_a in JRSSB paper

  rho_grid <- (1:1000)/1001

  alpha_grid <- alpha_scale * rho_grid / (1 - rho_grid )

  # logliks <- alpha_grid * mean_log_s +
  #   lgamma(alpha_grid) -
  #   p*lgamma(alpha_grid/p) +
  #   (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  #   # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)


  # logliks <- log(ddirichlet( t(matrix(s, p, 1000))  , t(matrix( rep(alpha_grid/p,p) , p , 1000)  ) ) ) +
  #   (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)

  logliks <- rep(NA, 1000)
  for(i in 1:1000){
    logliks[i] <- log(ddirichlet(s  , rep(alpha_grid[i]/p,p) ) ) +
      (alpha_a - 1)*log(rho_grid[i]) + (alpha_b-1)*log(1- rho_grid[i])
  }

  max_ll <- max(logliks)
  logsumexps <- max_ll + log(sum(exp( logliks  -  max_ll )))

  # print("logsumexps = ")
  # print(logsumexps)

  logliks <- exp(logliks - logsumexps)

  if(any(is.na(logliks))){
    print("logliks = ")
    print(logliks)

    print("logsumexps = ")
    print(logsumexps)

    print("mean_log_s = ")
    print(mean_log_s)

    print("lgamma(alpha_grid) = ")
    print(lgamma(alpha_grid))

    print("p*lgamma(alpha_grid/p) = ")
    print(p*lgamma(alpha_grid/p))

    print("(alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid) = ")
    print((alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid))

    print("max_ll = ")
    print(max_ll)

    print("s = ")
    print(s)


  }

  # print("logliks = ")
  # print(logliks)

  rho_ind <- sample.int(1000,size = 1, prob = logliks)


  return(alpha_grid[rho_ind])
}


update_sigma_mu <- function(trees, curr_sigmu2) {

  num_trees <- length(trees)
  mu_vec <- c()
  for(m in 1:length(trees)){

    mu_vec <- c(mu_vec,
                trees[[m]]$tree_matrix[, 'mu'])

  }

  mu_vec <- na.omit(mu_vec)

  # note Linero and Yang's sigma_mu corresponds to
  # num_trees times sigma_mu as defined in this package
  curr_s_mu <- sqrt(curr_sigmu2)

  prop_s_mu_minus2 <- rgamma(n = 1,
                             shape = length(mu_vec)/2,
                             rate = sum(mu_vec^2)/2)

  prop_s_mu <- sqrt(1/prop_s_mu_minus2)

  acceptprob <- (dcauchy(prop_s_mu, 0, 0.25/sqrt(num_trees))/dcauchy(curr_s_mu, 0, 0.25/sqrt(num_trees)))*
    (prop_s_mu/curr_s_mu)^3

  if(runif(1) < acceptprob){
    new_s_mu <- prop_s_mu
  }else{
    new_s_mu <- curr_s_mu
  }

  return(new_s_mu^2)
}


