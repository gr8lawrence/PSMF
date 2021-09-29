
# Update functions --------------------------------------------------------
## the following functions are for updating the parameters
update_C_TL <- function(C_old, M, G_old, beta) {

  # dimension
  n_cell_types = nrow(C_old)
  n_subjects = ncol(C_old)

  # constants
  num_matrix = t(G_old) %*% M + matrix(beta, n_cell_types, n_subjects)
  denom_matrix = t(G_old) %*% G_old %*% C_old + beta * matrix(rep(colSums(C_old), n_cell_types),
                                                              n_cell_types,
                                                              n_subjects,
                                                              byrow=TRUE)

  # updates (using the element-wise products)
  C_new = C_old * num_matrix/denom_matrix

  return(C_new)

}

update_G_TL <- function(G_old, M, C_new, G_0, n_markers, n_good_cell_types, alpha, xi) {

  # dimension
  n_genes = nrow(M)
  n_new_genes = n_genes - n_markers
  n_cell_types = nrow(C_new)
  n_bad_cell_types = n_cell_types - n_good_cell_types

  # constants
  Delta = rbind(
    matrix(
      rep(c(rep(1, n_good_cell_types), rep(0, n_bad_cell_types)), n_markers),
      n_markers,
      n_cell_types,
      byrow=TRUE
    ),
    matrix(
      0L,
      n_new_genes,
      n_cell_types
    )
  )
  Delta_c = 1 - Delta

  # update G
  G_new = matrix(0L, n_genes, n_cell_types)
  num_matrix= M %*% t(C_new) + alpha * Delta * G_0
  denom_matrix = G_old %*% C_new %*% t(C_new) + alpha * Delta * G_old  + xi * Delta_c * G_old
  G_new = G_old * num_matrix/denom_matrix

  return(G_new)

}

# Main algorithm ----------------------------------------------------------
## this section contains the main body of the algorithm
PSMF_solve <- function(M, G_0, G_init, C_init=NULL, n_markers, n_good_cell_types, alpha, xi, beta, max_iter=1e5, eps, method) {

  # parameters
  n_genes = nrow(M)
  n_subjects = ncol(M)
  n_cell_types = ncol(G_0)

  # initiate the counter
  n_iter = 0

  # inititate a residual vector
  res <- vector()

  # initiate matrices
  G_old = G_init
  if (method == 'auxiliary') {
    C_old = C_init
  } else {
    C_old = matrix(1e-8, n_cell_types, n_subjects)
  }

  # initiate the relative change in norm to infinity
  C_change = Inf
  G_change = Inf

  # start the algorithm
  while(n_iter < max_iter & min(C_change, G_change) >= eps) {

    if (n_iter %% 2 == 0) {
      # if n_iter is even, update C
      C_new = update_C(M=M,
                       G_old=G_old,
                       C_old=C_old,
                       beta=beta,
                       method=method)

      # get the relavitve change in norm for C
      C_change = norm(C_new - C_old, type="F")/norm(C_old, type="F")

      # replace the old C with new C
      C_old = C_new

      # calculate the residual
      res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")

    } else if (n_iter %% 2 == 1) {
      # if n_iter is odd, update G
      G_new = update_G(G_old=G_old,
                       M=M,
                       C_new=C_old,
                       G_0=G_0,
                       n_markers=n_markers,
                       n_good_cell_types=n_good_cell_types,
                       alpha=alpha,
                       xi=xi,
                       method=method)

      # get the relavitve change in norm for G
      G_change = norm(G_new - G_old, type="F")/norm(G_old, type="F")

      # replace the old G with new G
      G_old = G_new

      # calculate the residual
      res[n_iter + 1] = norm(M - G_old %*% C_old, type="F")

    }

    n_iter = n_iter + 1

  }

  #print(paste0("The algorithm terminated after finishing ", n_iter + 1, " iterations."))

  return(list(G_hat = G_old, C_hat = C_old, n_iter = n_iter, res_vec = res))

}
