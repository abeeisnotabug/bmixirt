generate_sinhsinh_stan_data <- function(N_levels, tol = sqrt(.Machine$double.eps), J) {
  m_t_max <-  log(2 * 2 / pi * log(2 * 2/pi * sqrt(.Machine$double.xmax)))

  max_num_absc <- ceiling(m_t_max / (2 / 2^(N_levels - 1)))

  m_abscissas <-
    m_weights <-
    m_weights_resc <- replicate(N_levels, rep(0, max_num_absc), simplify = FALSE)

  node_indices <- c(6L, 6L, 5L, 9L, 16L, 32L, 63L)

  for (i in 1:N_levels) {
    h <- 0.5^(i - 1)
    j <- 1
    arg <- h
    while (arg < m_t_max) {
      tmp <- pi/2 * sinh(arg)
      x <- sinh(tmp)
      m_abscissas[[i]][j] <- x
      m_weights[[i]][j] <- cosh(arg) * pi/2 * cosh(tmp)
      j <- j + 1
      if (i > 1) {
        arg <- arg + 2 * h
      } else {
        arg <- arg + h
      }
    }
    if (i > 2)
      m_weights_resc[[i]] <- m_weights[[i]] * h
    else
      m_weights_resc[[i]] <- m_weights[[i]]
  }
  abscissas_vec <- lapply(
    m_abscissas,
    function(level)
      t(
        sapply(
          level,
          rep,
          J
        )
      )
  )
  list(
    N_levels = N_levels,
    max_N_nodes = max_num_absc,
    abscissas = m_abscissas,
    abscissas_vec = abscissas_vec,
    weights = m_weights,
    weights_resc = m_weights_resc,
    node_indices = node_indices[seq_len(N_levels)],
    integrate = ifelse(tol < 1.0, 1L, 0L),
    tol = tol
  )
}

make_ucp_inits <- function(n_chains, C) {
  replicate(n_chains, list(uncond_class_prob = `dim<-`(rep(1/C, C), C)), simplify = FALSE)
}
