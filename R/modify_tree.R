# -------------------------------------------------------------------------#
# Description: this script contains functions that are used to perform the #
#              growing, pruning, changing, and swapping moves. It also has #
#              a function to initialise the trees to stumps                #
# -------------------------------------------------------------------------#

# 01. create_stump: initialises the trees to a stump
# 02. update_tree: calls the corresponding function associated to the move grow, prune, change, or swap.
# 03. grow_tree: grows a tree
# 04. prune_tree: prunes a tree
# 05. change_tree: changes the splitting rule that defines a pair of terminal nodes
# 06. swap_tree: exchanges the splitting rules that define two pair of terminal nodes

# Function to create all tree stumps ------------------------------------------------

create_stumps <- function(m,
                          y,
                          X) {
  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  #   Terminal (0 = no, 1 = yes)
  #   Child left
  #   Child right
  #   Node parents
  #   Split variable
  #   Split value
  #   mu
  #   Node size


  # Create holder for trees
  all_trees <- vector("list", length = m)
  # Initialisation of mu for the tree stumps
  mu_stump <- mean(y) #/ m
  # Loop through trees
  for (j in 1:m) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] <- vector("list", length = 2)
    # Give the elements names
    names(all_trees[[j]]) <- c("tree_matrix", "node_indices")
    # Create the two elements: first is a matrix and second is the assignment to node indices
    all_trees[[j]]$tree_matrix <- matrix(NA, ncol = 8, nrow = 1)
    all_trees[[j]]$node_indices <- rep(1, length(y))

    # Create column names
    colnames(all_trees[[j]][[1]]) <- c(
      "terminal",
      "child_left",
      "child_right",
      "parent",
      "split_variable",
      "split_value",
      "mu",
      "node_size"
    )

    # Set values for stump
    all_trees[[j]]$tree_matrix[1, ] <- c(1, NA, NA, NA, NA, NA, mu_stump, length(y))
  }
  return(all_trees)
}

# Function to update trees ------------------------------------------------

update_tree <- function(y, # Target variable
                        X, # Feature matrix
                        type = c(
                          "grow", # Grow existing tree
                          "prune", # Prune existing tree
                          "change" # Change existing tree - change split variable and value for an internal node
                        ),
                        curr_tree, # The current set of trees (not required if type is stump)
                        node_min_size, # The minimum size of a node to grow
                        s, # probability vector to be used during the growing process
                        max_bad_trees = 10) # Maximum of iterations to do to find a good tree with at least node_min_size observations in each leaf
{
  # Call the appropriate function to get the new tree
  new_tree <- switch(type,
                     grow = grow_tree(y, X, curr_tree, node_min_size, s, max_bad_trees),
                     prune = prune_tree(y, X, curr_tree, node_min_size),
                     change = change_tree(y, X, curr_tree, node_min_size, s, max_bad_trees)
  )

  return(new_tree)
}

# Grow_tree function ------------------------------------------------------

grow_tree <- function(y, X, curr_tree, node_min_size, s, max_bad_trees) {
  # Get the list of terminal nodes
  terminal_nodes <- as.numeric(which(curr_tree$tree_matrix[, "terminal"] == 1))

  # Find terminal node sizes
  terminal_node_size <- as.numeric(curr_tree$tree_matrix[terminal_nodes, "node_size"])

  # If a tree has no more leaves with at least 2 times node_min_size, then this tree cannot grow anymore.
  #   Then, return the current tree
  if (all(terminal_node_size < 2 * node_min_size)) {
    curr_tree$var <- 0
    return(curr_tree)
  }

  # Search for a good tree, which does exist using the case specific return above
  # However, only use a set maximum number of searches
  bad_tree <- TRUE
  n_bad_trees <- 0
  while(bad_tree){
    # Choose which node to split, set prob to zero for any nodes that are too small

    # Set up holder for new tree
    new_tree <- curr_tree

    # Add two extra rows to the tree in question
    new_tree$tree_matrix <- rbind(
      new_tree$tree_matrix,
      c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
      c(1, NA, NA, NA, NA, NA, NA, NA)
    )

    # node_to_split <- sample.vec(terminal_nodes, 1,
    #                             prob = as.integer(terminal_node_size >= 2 * node_min_size))
    if(length(terminal_nodes)==1){
      node_to_split <- terminal_nodes[1]
    }else{
      node_to_split <- sample(terminal_nodes, 1,
                                  prob = as.integer(terminal_node_size >= 2 * node_min_size))
    }



    # Choose a split variable using probability s (can be specified uniformly) from all columns
    # split_variable <- sample(1:ncol(X), 1, prob = s)
    split_variable <- sample.int(ncol(X), 1, prob = s)
    # available_values <- sort(unique(X[
    #   curr_tree$node_indices == node_to_split,
    #   split_variable
    # ]), na.last = TRUE)

    available_values <- collapse::funique(X[curr_tree$node_indices == node_to_split, split_variable], TRUE)

    # # If the number of unique values in the chosen node of the chosen covariate is less then 2 * node_min_size,
    # # then this is a bad tree choice!
    # n_values <- length(available_values)
    # if(n_values < 2 * node_min_size){
    #   n_bad_trees <- n_bad_trees + 1
    #   # If we reached the maximum of searches, stop searching for a good tree and return the current tree
    #   if(n_bad_trees >= max_bad_trees){
    #     curr_tree$var <- 0
    #     return(curr_tree)
    #   }
    # } else {
    #   bad_tree <- FALSE
    # }

    # can reduce the number of splits considered more by using
    # fnobs and checking if number above or below some available value is beyond the minimum node size

    if(length(available_values) == 1){
      new_split_value = available_values[1]
    } else if (length(available_values) == 2){
      new_split_value = available_values[2]
    }  else {
      # split_value = sample(available_values[-c(1,length(available_values))], 1)
      # split_value = resample(available_values[-c(1,length(available_values))])
      new_split_value = sample(available_values[-c(1)],1)
      # split_value = resample(available_values[-c(1,length(available_values))])
      # split_value = runif(1,available_values[2],available_values[length(available_values)])

    }


  # }

    # # Set up holder for new tree
    # new_tree <- curr_tree
    #
    # # Add two extra rows to the tree in question
    # new_tree$tree_matrix <- rbind(
    #   new_tree$tree_matrix,
    #   c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
    #   c(1, NA, NA, NA, NA, NA, NA, NA)
    # )

    # # Choose a split value such that there is enough observations to the left and right child
    # if (node_min_size == 0) {
    #   new_split_value <- resample(available_values)
    # } else if (node_min_size == 1) {
    #   new_split_value <- resample(available_values[-1])
    # } else {
    #   new_split_value <- resample(available_values[-c(1:node_min_size, (n_values - node_min_size + 2):n_values)])
    # }

    size_new_tree <- nrow(new_tree$tree_matrix)
    curr_parent <- new_tree$tree_matrix[node_to_split, "parent"] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split, 1:7] <- c(
      0, # No longer terminal
      size_new_tree - 1, # child_left is penultimate row
      size_new_tree, # child_right is last row
      curr_parent,
      split_variable,
      new_split_value,
      NA
    )

    #  Fill in the parents of these two nodes
    # new_tree$tree_matrix[size_new_tree,'parent'] = node_to_split
    new_tree$tree_matrix[c(size_new_tree - 1, size_new_tree), "parent"] <- node_to_split

    # Now call the fill function on this tree
    new_tree <- fill_tree_details(new_tree, X)

    # Store the covariate name to use it to update the Dirichlet prior of Linero (2016).
    new_tree$var <- split_variable

    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[,'node_size']) <= node_min_size)) {

      # print(" bad tree = ")
      # print(new_tree$tree_matrix)
      n_bad_trees = n_bad_trees + 1
    } else {
      bad_tree = FALSE
    }

    if(n_bad_trees >= max_bad_trees) {
      # print(" reached max_bad_trees = ")

      curr_tree$var = 0
      return(curr_tree)
    }
  }

  return(new_tree)
}

# Prune_tree function -----------------------------------------------------

prune_tree <- function(y, X, curr_tree, node_min_size) {
  # Create placeholder for new tree
  new_tree <- curr_tree

  if (nrow(new_tree$tree_matrix) == 1) { # No point in pruning a stump!
    new_tree$var <- 0
    return(new_tree)
  }

  # Get the list of terminal nodes
  terminal_nodes <- which(as.numeric(new_tree$tree_matrix[, "terminal"]) == 1)

  # Pick a random terminal node to prune, only pick nodes where both the left and right child are terminal!
  bad_node_to_prune <- TRUE # Start true
  while (bad_node_to_prune) {
    # Choose a random terminal node
    node_to_prune <- sample(terminal_nodes, 1)

    # Find the parent of this terminal node
    parent_pick <- as.numeric(new_tree$tree_matrix[node_to_prune, "parent"])
    var_pruned_nodes <- as.numeric(new_tree$tree_matrix[parent_pick, "split_variable"])

    # Get the two children of this parent
    child_left <- as.numeric(new_tree$tree_matrix[parent_pick, "child_left"])
    child_right <- as.numeric(new_tree$tree_matrix[parent_pick, "child_right"])

    # See whether either are terminal
    child_left_terminal <- as.numeric(new_tree$tree_matrix[child_left, "terminal"])
    child_right_terminal <- as.numeric(new_tree$tree_matrix[child_right, "terminal"])

    # If both are terminal then great
    if ((child_left_terminal == 1) & (child_right_terminal == 1)) {
      bad_node_to_prune <- FALSE # Have chosen a pair of terminal nodes so exist while loop
    } else {
      # Remove the node from the terminal node list to avoid sampling this node again
      terminal_nodes <- terminal_nodes[terminal_nodes != node_to_prune]
    }
  }

  # Delete these two rows from the tree matrix
  new_tree$tree_matrix <- new_tree$tree_matrix[-c(child_left, child_right), ,
                                               drop = FALSE
  ]
  # Make this node terminal again with no children or split values
  new_tree$tree_matrix[parent_pick, c(
    "terminal",
    "child_left",
    "child_right",
    "split_variable",
    "split_value"
  )] <- c(1, NA, NA, NA, NA)

  # If we're back to a stump no need to call fill_tree_details
  if (nrow(new_tree$tree_matrix) == 1) {
    new_tree$var <- var_pruned_nodes
    new_tree$node_indices <- rep(1, length(y))
  } else {
    # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
    if (node_to_prune <= nrow(new_tree$tree_matrix)) {
      # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
      bad_parents <- which(as.numeric(new_tree$tree_matrix[, "parent"]) >= node_to_prune)
      # Shift them back because of the two removed rows
      new_tree$tree_matrix[bad_parents, "parent"] <- as.numeric(new_tree$tree_matrix[bad_parents, "parent"]) - 2
      # Correct the parents and children for all the nodes after the node that was chosen to prune which have now shifted in the matrix
      for (j in node_to_prune:nrow(new_tree$tree_matrix)) {
        # Find the current parent
        curr_parent <- as.numeric(new_tree$tree_matrix[j, "parent"])
        # Find both the children of this node
        curr_children <- which(as.numeric(new_tree$tree_matrix[, "parent"]) == curr_parent)
        # Input these children back into the parent
        new_tree$tree_matrix[curr_parent, c("child_left", "child_right")] <- sort(curr_children)
      }
    }

    # Call the fill function on this tree
    new_tree <- fill_tree_details(new_tree, X)

    # Store the covariate name that was used in the splitting rule of the terminal nodes that were just pruned
    new_tree$var <- var_pruned_nodes
  }

  return(new_tree)
}

# change_tree function ----------------------------------------------------

change_tree <- function(y, X, curr_tree, node_min_size, s, max_bad_trees) {
  # Change a node means change out the split value and split variable of a second generation internal node (a node with only two children and no grandchildren etc.)

  # If current tree is a stump nothing to change
  if (nrow(curr_tree$tree_matrix) == 1) {
    curr_tree$var <- c(0, 0)
    return(curr_tree)
  }

  # Need to get the second generation internal nodes and the terminal nodes
  terminal_nodes <- which(as.numeric(curr_tree$tree_matrix[, "terminal"]) == 1)
  gen2_nodes <- get_gen2(curr_tree)

  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  count_bad_trees <- 0
  bad_trees <- TRUE
  while (bad_trees) {
    # Re-set the placeholder for the new tree
    new_tree <- curr_tree

    # Choose a second generation node to change. DO NOT USE STANDARD sample(...), this uses a different build-in feature for vectors of length 1!
    # node_to_change <- sample.vec(gen2_nodes, 1)

    if(length(gen2_nodes ==1)){
      node_to_change <- gen2_nodes[1]
    }else{
      node_to_change <- sample(gen2_nodes, 1)
    }

    # Get the covariate that will be changed
    var_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, "split_variable"])

    # Use the get_children function to get all the children (grandchildren, etc.) of this node
    all_children <- get_children(new_tree$tree_matrix, node_to_change)

    # Now find all the observations that fall in these children (these observations thus fall in the chosen node_to_change)
    use_node_indices <- !is.na(match(new_tree$node_indices, all_children))

    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree
    # new_split_variable <- sample(1:ncol(X), 1, prob = s)
    new_split_variable <- sample.int(ncol(X), 1, prob = s)
    # available_values <- sort(unique(X[
    #   use_node_indices,
    #   new_split_variable
    # ]))
    available_values <- collapse::funique(X[use_node_indices,new_split_variable], TRUE)

    n_values <- length(available_values)

    # # If there are not enough available values to at least assign node_min_size amount of observations to the left
    # #   and right, then it is already a bad tree
    # if (n_values < 2 * node_min_size) {
    #   count_bad_trees <- count_bad_trees + 1
    #   if (count_bad_trees >= max_bad_trees) {
    #     curr_tree$var <- c(0, 0)
    #     return(curr_tree)
    #   }
    #   # If this a bad tree, skip the rest of this while iteration
    #   next
    # } else {
    #   # Prevent a bad tree to remove the first and (last-1) node_min_size values
    #   if (node_min_size == 0) {
    #     new_split_value <- resample(available_values)
    #   } else if (node_min_size == 1) {
    #     new_split_value <- resample(available_values[-1])
    #   } else {
    #     # Does not work if node_min_size is 0 or 1, hence the ifs
    #     new_split_value <- resample(available_values[-c(1:node_min_size, (n_values - node_min_size + 2):n_values)])
    #   }
    # }

    if (length(available_values) == 1){
      new_split_value = available_values[1]
      new_tree$var = c(var_changed_node, new_split_variable)
    } else if (length(available_values) == 2){
      new_split_value = available_values[2]
      new_tree$var = c(var_changed_node, new_split_variable)
    } else {
      # new_split_value = resample(available_values[-c(1,length(available_values))], 1)
      # new_split_value = sample(available_values[-c(1,length(available_values))])
      new_split_value = sample(available_values[-c(1)],1)
      # new_split_value = runif(1,available_values[2],available_values[length(available_values)])

    }

    # Update the tree details
    new_tree$tree_matrix[
      node_to_change,
      c(
        "split_variable",
        "split_value"
      )
    ] <- c(
      new_split_variable,
      new_split_value
    )

    # Update the tree node indices
    new_tree <- fill_tree_details(new_tree, X)

    # Store the covariate name that was used in the splitting rule of the terminal node that was just changed
    new_tree$var <- c(var_changed_node, new_split_variable)

    # Another check for bad tree since the children of the changed node can be internal nodes with now fewer then
    #   node_min_size unique observations
    if (any(as.numeric(new_tree$tree_matrix[terminal_nodes, "node_size"]) < node_min_size)) {
      count_bad_trees <- count_bad_trees + 1
      if (count_bad_trees == max_bad_trees) {
        curr_tree$var <- c(0, 0)
        return(curr_tree)
      }
    } else {
      bad_trees <- FALSE
    }
  }

  # # this if statement is probably unnecessary, unless "next" occurs in last while looop iteration
  # if(bad_trees == TRUE){
  #   curr_tree$var <- c(0, 0)
  # }

  return(new_tree)
}

change_tree_global <- function(y, X, curr_tree, node_min_size, s, max_bad_trees) {
  # Change a node means change out the split value and split variable of any internal node. Need to make sure that this does not produce a bad tree (i.e. zero terminal nodes)

  # If current tree is a stump nothing to change
  if (nrow(curr_tree$tree_matrix) == 1) {
    curr_tree$var <- c(0, 0)
    return(curr_tree)
  }

  # Create a holder for the new tree
  new_tree <- curr_tree

  # Need to get the internal nodes
  internal_nodes <- which(as.numeric(new_tree$tree_matrix[, "terminal"]) == 0)
  terminal_nodes <- which(as.numeric(new_tree$tree_matrix[, "terminal"]) == 1)

  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  count_bad_trees <- 0
  bad_trees <- TRUE
  while (bad_trees) {
    # Re-set the tree
    new_tree <- curr_tree

    # choose an internal node to change
    node_to_change <- sample(internal_nodes, 1)

    # Get the covariate that will be changed
    var_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, "split_variable"])

    # Use the get_children function to get all the children (grandchildren, etc.) of this node
    all_children <- get_children(new_tree$tree_matrix, node_to_change)

    # Now find all the observations that fall in these children (these observations thus fall in the chosen node_to_change)
    use_node_indices <- !is.na(match(new_tree$node_indices, all_children))

    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree
    # new_split_variable <- sample(1:ncol(X), 1, prob = s)
    new_split_variable <- sample.int(ncol(X), 1, prob = s)
    # available_values <- sort(unique(X[
    #   use_node_indices,
    #   new_split_variable
    # ]))
    available_values <- collapse::funique(X[use_node_indices,new_split_variable], TRUE)

    n_values <- length(available_values)

    # If there are not enough available values to at least assign node_min_size amount of observations to the left
    #   and right, then it is already a bad tree
    if (n_values < 2 * node_min_size) {
      count_bad_trees <- count_bad_trees + 1
      if (count_bad_trees == max_bad_trees) {
        curr_tree$var <- c(0, 0)
        return(curr_tree)
      }
      next
    } else {
      # Prevent a bad tree to remove the first and (last-1) node_min_size values
      if (node_min_size == 0) {
        new_split_value <- resample(available_values)
      } else if (node_min_size == 1) {
        new_split_value <- resample(available_values[-1])
      } else {
        # Does not work if node_min_size is 0 or 1, hence the ifs
        new_split_value <- resample(available_values[-c(1:node_min_size, (n_values - node_min_size + 2):n_values)])
      }
    }

    # Update the tree details
    new_tree$tree_matrix[
      node_to_change,
      c(
        "split_variable",
        "split_value"
      )
    ] <- c(
      new_split_variable,
      new_split_value
    )

    # Update the tree node indices
    new_tree <- fill_tree_details(new_tree, X)

    # Store the covariate name that was used in the splitting rule of the terminal node that was just changed
    new_tree$var <- c(var_changed_node, new_split_variable)

    # Another check for bad tree since the children of the changed node can be internal nodes with now fewer then
    #   node_min_size observations
    if (any(as.numeric(new_tree$tree_matrix[terminal_nodes, "node_size"]) < node_min_size)) {
      count_bad_trees <- count_bad_trees + 1
    } else {
      bad_trees <- FALSE
    }
    if (count_bad_trees == max_bad_trees) {
      curr_tree$var <- c(0, 0)
      return(curr_tree)
    }
  }

  return(new_tree)
}
