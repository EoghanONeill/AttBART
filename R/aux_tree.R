# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. get_number_distinct_cov: given a tree, it returns the number of distinct covariates used to create its structure

# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices
  node_indices = rep(1, nrow(X))

  # print("nrow(tree_matrix) = ")
  # print(nrow(tree_matrix))


  if(nrow(tree_matrix) >1){
    # For all but the top row, find the number of observations falling into each one
    for(i in 2:nrow(tree_matrix)) {

      # Get the parent
      curr_parent = as.numeric(tree_matrix[i,'parent'])

      # Find the split variable and value of the parent
      split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
      split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

      # Find whether it's a left or right terminal node
      left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                             'left', 'right')
      if(left_or_right == 'left') {
        # If left use less than condition
        new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
        node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
      } else {
        # If right use greater than condition
        new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
        node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
      }
    } # End of loop through table
  }


  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]


  # THIS FUNCTION MUST BE COMPLEtely re-WRITTEN FOR AttBART
  # Attention weights must be calculated and normalized

  # Need Separate function for test data with additional input corresponding to test covariate matrix


  # LATER INCLUDE EPSILON AND w as INPUTS

  # TO calculate Weights
    # Find Leaf for each training observation
    # (Find leaf for test observation(s) if applicable)
    # Obtain leaf covariate means (vectors)
    # For each (training or test) observation, subtract the covariate mean ofthe relevant leaf
    # Square the difference, then divide by 2, then take the exponential
    # normalize to sum to 1 (softmax)
    # later will account for epsilon and w
    # set equal to attention weights

  # For each tree, multiply each individuals prediciton by relevant weight

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # If one tree, then attention weight must be 1, so the predicitons are unaffected

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


    if(any(is.na(predictions))){

      print("trees = ")
      print(trees)

      print("predictions = ")
      print(predictions)


      stop("Error in get_predictions. NA prediction values")
    }


  } else {
    # # Do a recursive call to the function
    # partial_trees = trees
    # partial_trees[[1]] = NULL # Blank out that element of the list
    # predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
    #   get_predictions(partial_trees, X,
    #                   single_tree = length(partial_trees) == 1)
    # #single_tree = !is.null(names(partial_trees)))
    # # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)


    # cannot proceed with recursion, instead must save weights and leaf means

    # print("length(trees) = ")
    # print(length(trees))
    #
    # print("nrow(X) = ")
    # print(nrow(X))

    # matrix of tree predictions
    treepreds <- matrix(NA,
                        nrow = nrow(X),
                        ncol = length(trees))

    #matrix of observation and tree-specific unnormalized weights
    expdistweights <- matrix(NA,
                        nrow = nrow(X),
                        ncol = length(trees))


    for(tree_ind in 1:length(trees)){

      # save unweighted predictions for tree
      treepreds[, tree_ind] <- get_predictions(trees[[tree_ind]], X, single_tree = TRUE)

      temptree <- trees[[tree_ind]]
      unique_node_indices = unique(temptree$node_indices)
      # Get the node indices for the current X matrix
      # print("temptree = ")
      # print(temptree)
      # print("X = ")
      # print(X)


      curr_X_node_indices = fill_tree_details(temptree, X)$node_indices

      # create mean vectors for each node
      # let each column correspond to a node
      # leafmeanmat <- matrix(NA, nrow = ncol(X), ncol = length(unique_node_indices))

      # for each X row, subtract the relevant leaf mean
      diffXmat <- X

      for(i in 1:length(unique_node_indices)) {
        #column (variable) means of subset of X matrix containing observations in leaf i
        # leafmeanmat[,i] = apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )
        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   leafmeanmat[,i]

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )


        diffXmat[curr_X_node_indices == unique_node_indices[i],] =
          sweep(diffXmat[curr_X_node_indices == unique_node_indices[i],, drop = FALSE],
                2,
                apply(X[curr_X_node_indices == unique_node_indices[i],, drop = FALSE],2,mean ))

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        # apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )

      }

      # square differences (and arbitrarily divide by 2)
      diffXmat <- (diffXmat^2)
      # obtain row sums of squared differences
      distvec <- apply(diffXmat,1,sum)

      # arbitrarily divide by 2 and take the exponential
      distvec <- exp( -1*distvec/2 )


      #

      # Now for each X row, subtract the relevant leaf mean
      # diffXmat <- X
      #
      # # Now loop through all node indices to fill in details
      # for(i in 1:length(unique_node_indices)) {
      #   diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
      #     leafmeanmat[,i]
      # }


      expdistweights[, tree_ind] <- distvec # ENTER Weight vector here

    }

    #normalize distances to obtain softmax weights
    # divide each row by row mean
    smaxweights <- t(apply(expdistweights,1, function(x) x/sum(x)))

    # LATER ACCOUNT FOR EPSILON AND W HERE
    # for now, just use smax weights as attention weights
    #
    #
    #
    #
    #

    Aweights <- smaxweights

    # multiply predictions and weights element-wise
    # then take the row sums
    predictions <-  apply(treepreds*smaxweights,1,sum)

    if(any(is.na(predictions))){

      print("treepreds = ")
      print(treepreds)

      print("smaxweights = ")
      print(smaxweights)


      stop("Error in get_predictions. NA prediction values")
    }

  }

  if(any(is.na(predictions))){

    print("treepreds = ")
    print(treepreds)

    print("smaxweights = ")
    print(smaxweights)


    stop("Error in get_predictions. NA prediction values")
  }
  return(predictions)
}



# Get predictions test ---------------------------------------------------------

get_predictions_test = function(trees, X, single_tree = FALSE, xtrain) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]


  # THIS FUNCTION MUST BE COMPLEtely re-WRITTEN FOR AttBART
  # Attention weights must be calculated and normalized

  # Need Separate function for test data with additional input corresponding to test covariate matrix


  # LATER INCLUDE EPSILON AND w as INPUTS

  # TO calculate Weights
  # Find Leaf for each training observation
  # (Find leaf for test observation(s) if applicable)
  # Obtain leaf covariate means (vectors)
  # For each (training or test) observation, subtract the covariate mean ofthe relevant leaf
  # Square the difference, then divide by 2, then take the exponential
  # normalize to sum to 1 (softmax)
  # later will account for epsilon and w
  # set equal to attention weights

  # For each tree, multiply each individuals prediciton by relevant weight

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # If one tree, then attention weight must be 1, so the predicitons are unaffected

    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      # unique_node_indices = unique(trees$node_indices)

      temptreemat <- trees$tree_matrix
      unique_node_indices <- which(temptreemat[,1] ==1)

      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {

        # print("line 322. curr_X_node_indices = ")
        # print(curr_X_node_indices)
        #
        # print("unique_node_indices = ")
        # print(unique_node_indices)
        #
        # print("temptreemat = ")
        # print(temptreemat)
        #
        # print("i = ")
        # print(i)


        if(sum(curr_X_node_indices == unique_node_indices[i])==0){
          #next #no test observations in leaf
        }else{
          predictions[curr_X_node_indices == unique_node_indices[i]] =
            trees$tree_matrix[unique_node_indices[i], 'mu']
        }
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees


    if(any(is.na(predictions))){

      print("trees = ")
      print(trees)

      print("predictions = ")
      print(predictions)


      stop("Error in get_predictions. NA prediction values")
    }


  } else {
    # # Do a recursive call to the function
    # partial_trees = trees
    # partial_trees[[1]] = NULL # Blank out that element of the list
    # predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
    #   get_predictions(partial_trees, X,
    #                   single_tree = length(partial_trees) == 1)
    # #single_tree = !is.null(names(partial_trees)))
    # # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)


    # cannot proceed with recursion, instead must save weights and leaf means

    # print("length(trees) = ")
    # print(length(trees))
    #
    # print("nrow(X) = ")
    # print(nrow(X))

    # matrix of tree predictions
    treepreds <- matrix(NA,
                        nrow = nrow(X),
                        ncol = length(trees))

    #matrix of observation and tree-specific unnormalized weights
    expdistweights <- matrix(NA,
                             nrow = nrow(X),
                             ncol = length(trees))


    for(tree_ind in 1:length(trees)){

      # save unweighted predictions for tree
      treepreds[, tree_ind] <- get_predictions(trees[[tree_ind]], X, single_tree = TRUE)

      temptree <- trees[[tree_ind]]
      # unique_node_indices = unique(temptree$node_indices)

      temptreemat <- temptree$tree_matrix
      unique_node_indices <- which(temptreemat[,1] ==1)

      # Get the node indices for the current X matrix
      # print("temptree = ")
      # print(temptree)
      # print("X = ")
      # print(X)


      curr_X_node_indices = fill_tree_details(temptree, X)$node_indices
      curr_X_node_indices_train = fill_tree_details(temptree, xtrain)$node_indices

      # create mean vectors for each node
      # let each column correspond to a node
      # leafmeanmat <- matrix(NA, nrow = ncol(X), ncol = length(unique_node_indices))

      # for each X row, subtract the relevant leaf mean
      diffXmat <- matrix(NA,
                         nrow = nrow(X),
                         ncol = ncol(X))

      for(i in 1:length(unique_node_indices)) {
        #column (variable) means of subset of X matrix containing observations in leaf i
        # leafmeanmat[,i] = apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )
        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   leafmeanmat[,i]

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )

        # print("line 428. curr_X_node_indices = ")
        # print(curr_X_node_indices)
        #
        # print("unique_node_indices = ")
        # print(unique_node_indices)
        #
        # print("i = ")
        # print(i)
        #
        # print("temptreemat = ")
        # print(temptreemat)


        if(sum(curr_X_node_indices == unique_node_indices[i])==0){
          # DO NOTHING
          # no test observations in relevant leaf
        }else{

          # print("curr_X_node_indices = ")
          # print(curr_X_node_indices)
          # print("curr_X_node_indices_train = ")
          # print(curr_X_node_indices_train)
          #
          # print("unique_node_indices[i] = ")
          # print(unique_node_indices[i])
          #
          #
          # print("diffXmat[curr_X_node_indices == unique_node_indices[i],  ]= ")
          # print(diffXmat[curr_X_node_indices == unique_node_indices[i],  ])

          diffXmat[curr_X_node_indices == unique_node_indices[i],  ] =
            sweep(X[curr_X_node_indices == unique_node_indices[i], , drop = FALSE],
                  2,
                  apply(xtrain[curr_X_node_indices_train == unique_node_indices[i],, drop = FALSE],2,mean ))
        }
        # diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        # apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )

      }

      # square differences (and arbitrarily divide by 2)
      diffXmat <- (diffXmat^2)
      # obtain row sums of squared differences
      distvec <- apply(diffXmat,1,sum)

      # arbitrarily divide by 2 and take the exponential
      distvec <- exp( -1*distvec/2 )


      #

      # Now for each X row, subtract the relevant leaf mean
      # diffXmat <- X
      #
      # # Now loop through all node indices to fill in details
      # for(i in 1:length(unique_node_indices)) {
      #   diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
      #     leafmeanmat[,i]
      # }


      expdistweights[, tree_ind] <- distvec # ENTER Weight vector here

    }

    #normalize distances to obtain softmax weights
    # divide each row by row mean
    smaxweights <- t(apply(expdistweights,1, function(x) x/sum(x)))

    # LATER ACCOUNT FOR EPSILON AND W HERE
    # for now, just use smax weights as attention weights
    #
    #
    #
    #
    #

    Aweights <- smaxweights

    # multiply predictions and weights element-wise
    # then take the row sums
    predictions <-  apply(treepreds*smaxweights,1,sum)

    if(any(is.na(predictions))){

      print("treepreds = ")
      print(treepreds)

      print("smaxweights = ")
      print(smaxweights)


      stop("Error in get_predictions. NA prediction values")
    }

  }

  if(any(is.na(predictions))){

    print("treepreds = ")
    print(treepreds)

    print("smaxweights = ")
    print(smaxweights)


    stop("Error in get_predictions. NA prediction values")
  }
  return(predictions)
}






# Get predictions without tree ---------------------------------------------------------

get_predictions_drop = function(trees, X, single_tree = FALSE, drop_ind) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]


  # THIS FUNCTION MUST BE COMPLEtely re-WRITTEN FOR AttBART
  # Attention weights must be calculated and normalized

  # Need Separate function for test data with additional input corresponding to test covariate matrix


  # LATER INCLUDE EPSILON AND w as INPUTS

  # TO calculate Weights
  # Find Leaf for each training observation
  # (Find leaf for test observation(s) if applicable)
  # Obtain leaf covariate means (vectors)
  # For each (training or test) observation, subtract the covariate mean ofthe relevant leaf
  # Square the difference, then divide by 2, then take the exponential
  # normalize to sum to 1 (softmax)
  # later will account for epsilon and w
  # set equal to attention weights

  # For each tree, multiply each individuals prediciton by relevant weight

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # If one tree, then attention weight must be 1, so the predicitons are unaffected

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


    if(any(is.na(predictions))){

      print("trees = ")
      print(trees)

      print("predictions = ")
      print(predictions)


      stop("Error in get_predictions. NA prediction values")
    }


  } else {
    # # Do a recursive call to the function
    # partial_trees = trees
    # partial_trees[[1]] = NULL # Blank out that element of the list
    # predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
    #   get_predictions(partial_trees, X,
    #                   single_tree = length(partial_trees) == 1)
    # #single_tree = !is.null(names(partial_trees)))
    # # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)


    # cannot proceed with recursion, instead must save weights and leaf means

    # matrix of tree predictions
    treepreds <- matrix(NA,
                        nrow = nrow(X),
                        ncol = length(trees))

    #matrix of observation and tree-specific unnormalized weights
    expdistweights <- matrix(NA,
                             nrow = nrow(X),
                             ncol = length(trees))


    for(tree_ind in 1:length(trees)){

      # print("tree_ind =  ")
      # print(tree_ind)
      #
      # print("drop_ind =  ")
      # print(drop_ind)

      if(tree_ind != drop_ind){
        # save unweighted predictions for tree
        treepreds[, tree_ind] <- get_predictions(trees[[tree_ind]], X, single_tree = TRUE)
      }

      temptree <- trees[[tree_ind]]
      unique_node_indices = unique(temptree$node_indices)
      # Get the node indices for the current X matrix
      # print("temptree = ")
      # print(temptree)
      # print("X = ")
      # print(X)


      curr_X_node_indices = fill_tree_details(temptree, X)$node_indices

      # create mean vectors for each node
      # let each column correspond to a node
      # leafmeanmat <- matrix(NA, nrow = ncol(X), ncol = length(unique_node_indices))

      # for each X row, subtract the relevant leaf mean
      diffXmat <- X

      for(i in 1:length(unique_node_indices)) {
        #column (variable) means of subset of X matrix containing observations in leaf i
        # leafmeanmat[,i] = apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )
        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   leafmeanmat[,i]

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )


        diffXmat[curr_X_node_indices == unique_node_indices[i],] =
          sweep(diffXmat[curr_X_node_indices == unique_node_indices[i],, drop = FALSE],
                2,
                apply(X[curr_X_node_indices == unique_node_indices[i],, drop = FALSE],2,mean ))

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        # apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )

      }

      # square differences (and arbitrarily divide by 2)
      diffXmat <- (diffXmat^2)
      # obtain row sums of squared differences
      distvec <- apply(diffXmat,1,sum)

      # arbitrarily divide by 2 and take the exponential
      distvec <- exp( -1*distvec/2 )


      #

      # Now for each X row, subtract the relevant leaf mean
      # diffXmat <- X
      #
      # # Now loop through all node indices to fill in details
      # for(i in 1:length(unique_node_indices)) {
      #   diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
      #     leafmeanmat[,i]
      # }


      expdistweights[, tree_ind] <- distvec # ENTER Weight vector here

    }

    #normalize distances to obtain softmax weights
    # divide each row by row mean
    smaxweights <- t(apply(expdistweights,1, function(x) x/sum(x)))

    # LATER ACCOUNT FOR EPSILON AND W HERE
    # for now, just use smax weights as attention weights
    #
    #
    #
    #
    #

    Aweights <- smaxweights

    # multiply predictions and weights element-wise
    # then take the row sums
    predictions <-  apply(treepreds[,-drop_ind]*smaxweights[,-drop_ind],1,sum)

    if(any(is.na(predictions))){

      print("treepreds = ")
      print(treepreds)

      print("smaxweights = ")
      print(smaxweights)


      stop("Error in get_predictions. NA prediction values")
    }

  }

  if(any(is.na(predictions))){

    print("treepreds = ")
    print(treepreds)

    print("smaxweights = ")
    print(smaxweights)


    stop("Error in get_predictions. NA prediction values")
  }
  return(predictions)
}





# Get attention ---------------------------------------------------------

get_attention = function(trees, X) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]


  # THIS FUNCTION MUST BE COMPLEtely re-WRITTEN FOR AttBART
  # Attention weights must be calculated and normalized

  # Need Separate function for test data with additional input corresponding to test covariate matrix


  # LATER INCLUDE EPSILON AND w as INPUTS

  # TO calculate Weights
  # Find Leaf for each training observation
  # (Find leaf for test observation(s) if applicable)
  # Obtain leaf covariate means (vectors)
  # For each (training or test) observation, subtract the covariate mean ofthe relevant leaf
  # Square the difference, then divide by 2, then take the exponential
  # normalize to sum to 1 (softmax)
  # later will account for epsilon and w
  # set equal to attention weights

  # For each tree, multiply each individuals prediciton by relevant weight

  # Normally trees will be a list of lists but just in case
  # if(single_tree) {
  #
  #   # If one tree, then attention weight must be 1, so the predicitons are unaffected
  #
  #   # Deal with just a single tree
  #   if(nrow(trees$tree_matrix) == 1) {
  #     predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
  #   } else {
  #     # Loop through the node indices to get predictions
  #     predictions = rep(NA, nrow(X))
  #     unique_node_indices = unique(trees$node_indices)
  #     # Get the node indices for the current X matrix
  #     curr_X_node_indices = fill_tree_details(trees, X)$node_indices
  #     # Now loop through all node indices to fill in details
  #     for(i in 1:length(unique_node_indices)) {
  #       predictions[curr_X_node_indices == unique_node_indices[i]] =
  #         trees$tree_matrix[unique_node_indices[i], 'mu']
  #     }
  #   }
  #   # More here to deal with more complicated trees - i.e. multiple trees
  # } else {
    # Do a recursive call to the function
    # partial_trees = trees
    # partial_trees[[1]] = NULL # Blank out that element of the list
    # predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
    #   get_predictions(partial_trees, X,
    #                   single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)


    # cannot proceed with recursion, instead must save weights and leaf means

    # matrix of tree predictions
    # treepreds <- matrix(NA,
    #                     nrow = nrow(X),
    #                     ncol = length(trees))

    #matrix of observation and tree-specific unnormalized weights
    expdistweights <- matrix(NA,
                             nrow = nrow(X),
                             ncol = length(trees))


    for(tree_ind in 1:length(trees)){

      # save unweighted predictions for tree
      # treepreds[, tree_ind] <- get_predictions(trees[[tree_ind]], X, single_tree = TRUE)

      temptree <- trees[[tree_ind]]
      unique_node_indices = unique(temptree$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(temptree, X)$node_indices

      # create mean vectors for each node
      # let each column correspond to a node
      # leafmeanmat <- matrix(NA, nrow = ncol(X), ncol = length(unique_node_indices))

      # for each X row, subtract the relevant leaf mean
      diffXmat <- X

      for(i in 1:length(unique_node_indices)) {
        #column (variable) means of subset of X matrix containing observations in leaf i
        # leafmeanmat[,i] = apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )
        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   leafmeanmat[,i]

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        #   apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )


        diffXmat[curr_X_node_indices == unique_node_indices[i],] =
          sweep(diffXmat[curr_X_node_indices == unique_node_indices[i],],
                2,
                apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean ))

        # diffXmat[curr_X_node_indices == unique_node_indices[i],] -
        # apply(X[curr_X_node_indices == unique_node_indices[i],],2,mean )

      }

      # square differences (and arbitrarily divide by 2)
      diffXmat <- (diffXmat^2)
      # obtain row sums of squared differences
      distvec <- apply(diffXmat,1,sum)

      # arbitrarily divide by 2 and take the exponential
      distvec <- exp( -1*distvec/2 )


      #

      # Now for each X row, subtract the relevant leaf mean
      # diffXmat <- X
      #
      # # Now loop through all node indices to fill in details
      # for(i in 1:length(unique_node_indices)) {
      #   diffXmat[curr_X_node_indices == unique_node_indices[i],] = diffXmat[curr_X_node_indices == unique_node_indices[i],] -
      #     leafmeanmat[,i]
      # }


      expdistweights[, tree_ind] <- distvec # ENTER Weight vector here

    }

    #normalize distances to obtain softmax weights
    # divide each row by row sum
    smaxweights <- t(apply(expdistweights,1, function(x) x/sum(x)))

    # LATER ACCOUNT FOR EPSILON AND W HERE
    # for now, just use smax weights as attention weights
    #
    #
    #
    #
    #

    Aweights <- smaxweights

    # multiply predictions and weights element-wise
    # then take the row sums
    # predictions <-  apply(treepreds*smaxweights,1,sum)



  # }

  return(Aweights)
}







# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

update_s = function(var_count, p, alpha_s){
  s_ = rdirichlet(1, alpha_s/p + var_count)
  return(s_)
}

get_number_distinct_cov <- function(tree){

  # Select the rows that correspond to internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
  # Get the covariates that are used to define the splitting rules
  num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))

  return(num_distinct_cov)
}
