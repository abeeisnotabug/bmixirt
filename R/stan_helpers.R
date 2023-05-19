generate_sinhsinh_stan_data <- function(N_levels, tol = sqrt(.Machine$double.eps), J = FALSE) {
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
  if (!isFALSE(J))
    abscissas_vec <- lapply(m_abscissas,
                            function(level) t(sapply(level, rep, J)))

  list(
    N_levels = N_levels,
    max_N_nodes = max_num_absc,
    abscissas = m_abscissas,
    if (!isFALSE(J)) abscissas_vec = abscissas_vec,
    weights = m_weights,
    weights_resc = m_weights_resc,
    node_indices = node_indices[seq_len(N_levels)],
    integrate = ifelse(tol < 1.0, 1L, 0L),
    tol = tol
  )
}

mccirt_r_itsl_v0_initfun <- function(K, J, C, N_R, N_T, N_I)
  list(uncond_classprob = rep(1/C, C),
       threshold_c = abind::abind(lapply(1:C, function(bnd) array(rep(seq(-bnd, bnd, length.out = J), each = K), dim = c(1, K, J))), along = 1),
       log_loading_I_c = rep(0, K),
       log_loading_R_c = rep(0, K),
       log_loading_T_c = rep(0, K),
       intercept = rep(0, C),
       log_scaling_I = 0,
       log_scaling_R = rep(0, C),
       log_scaling_T = 0,
       ability_I = rep(0, N_I),
       ability_R = rep(0, N_R),
       ability_T = rep(0, N_T))

compile_model <- function(cmdstan_path = cmdstanr::cmdstan_path(), path_prefix = file.path("data", "sim", "li53vet"), path, wait = FALSE)
  system(sprintf("cd %s\nmake %s",
                 cmdstan_path, file.path(path_prefix, path)),
         wait = wait)
