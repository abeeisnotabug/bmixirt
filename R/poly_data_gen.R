dmnorm <- function(x, ucp, means, sds) {
  dnorm(sapply(seq_along(means), function(i) (x - means[i]) / sds[i])) %*% (ucp / sds)
}

expected_entropy <- function(ucp, means, sds, order = 5) {
  this_log_dmnorm <- function(x) {
    log(dmnorm(x, ucp = ucp, means = means, sds = sds))
  }

  C <- length(ucp)

  first_summands <- sapply(
    1:C,
    function(c)
      spatstat.core::gauss.hermite(this_log_dmnorm, mu = means[c], sd = sds[c], order = order)
  )

  second_summands <- 0.5 + log(sqrt(2 * pi)) + log(sds) - log(ucp)

  1 - t(first_summands + second_summands) %*% ucp / log(C)
}

grm_cum_prob <- function(ability, alpha_loading, betas_thresholds, irt) {
  if (irt)
    t(1 / (1 + exp(alpha_loading * outer(betas_thresholds, ability, "-"))))
  else
    t(1 / (1 + exp(outer(betas_thresholds, alpha_loading * ability, "-"))))
}

grm_cat_prob <- function(ability, alpha_loading, betas_thresholds, irt) {
  cat_probs <- cbind(
    rep(1, length(ability)),
    grm_cum_prob(ability, alpha_loading, betas_thresholds, irt)
  )

  cbind(
    -t(diff(t(cat_probs))),
    cat_probs[, ncol(cat_probs)]
  )
}

gpcm_cat_prob <- function(ability, alpha_loading, betas_thresholds, irt, normalized = TRUE) {
  num <- cbind(
    rep(1, length(ability)),
    if (irt)
      exp(alpha_loading * sweep(ability %o% seq_along(betas_thresholds), 2, cumsum(betas_thresholds), "-"))
    else
      exp(sweep(alpha_loading * ability %o% seq_along(betas_thresholds), 2, cumsum(betas_thresholds), "-"))
  )

  if (normalized)
    sweep(num, 1, rowSums(num), "/")
  else
    num
}

gpcm_cum_prob <- function(ability, alpha_loading, betas_thresholds, irt) {
  cat_probs <- gpcm_cat_prob(ability, alpha_loading, betas_thresholds, irt, normalized = TRUE)

  t(apply(cat_probs[, ncol(cat_probs):2, drop = FALSE], 1, cumsum))[, (ncol(cat_probs) - 1):1, drop = FALSE]
}

gen_pirt <- function(ability, alpha_loading, betas_thresholds, irt, model) {
  if (model == "gpcm") {
    probs <- gpcm_cat_prob(ability, alpha_loading, betas_thresholds, irt, normalized = FALSE)
  } else if (model == "grm") {
    probs <- grm_cat_prob(ability, alpha_loading, betas_thresholds, irt)
  }

  extraDistr::rcat(length(ability), probs)
}

convert_pars_to_other_model <- function(
  alpha_loading,
  betas_thresholds,
  irt,
  from,
  feature,
  rtol = 1e-12,
  atol = 1e-16,
  ctol = 1e-16
) {
  grm_diff_fun <- function(alpha_loading, given_betas_thresholds, irt, mode = "find_parameters") {
    function(betas_thresholds_to_find) {
      cum_prob <- if (mode == "find_parameters") {
        cbind(
          rep(1, length(given_betas_thresholds)),
          grm_cum_prob(given_betas_thresholds, alpha_loading, betas_thresholds_to_find, irt)
        )
      } else if (mode == "find_intersection") {
        cbind(
          rep(1, length(given_betas_thresholds)),
          grm_cum_prob(betas_thresholds_to_find, alpha_loading, given_betas_thresholds, irt)
        )
      }

      cat_prob <- cbind(
        -t(diff(t(cum_prob))),
        cum_prob[, ncol(cum_prob)]
      )

      diffs <- -t(diff(t(cat_prob)))
      diag(diffs)
    }
  }

  gpcm_half_fun <- function(alpha_loading, given_betas_thresholds, irt, mode = "find_parameters") {
    function(betas_thresholds_to_find) {
      if (mode == "find_parameters") {
        diag(gpcm_cum_prob(given_betas_thresholds, alpha_loading, betas_thresholds_to_find, irt)) - 0.5
      } else if (mode == "find_intersection") {
        diag(gpcm_cum_prob(betas_thresholds_to_find, alpha_loading, given_betas_thresholds, irt)) - 0.5
      }
    }
  }

  if (from == "gpcm") {
    if (feature == "half") {
      rootSolve::multiroot(
        gpcm_half_fun(alpha_loading, betas_thresholds, irt, mode = "find_intersection"),
        betas_thresholds,
        rtol = rtol,
        atol = atol,
        ctol = ctol
      )$root
    } else if (feature == "intersection") {
      rootSolve::multiroot(
        grm_diff_fun(alpha_loading, betas_thresholds, irt, mode = "find_parameters"),
        betas_thresholds,
        rtol = rtol,
        atol = atol,
        ctol = ctol
      )$root
    }
  } else if (from == "grm") {
    if (feature == "intersection") {
      rootSolve::multiroot(
        grm_diff_fun(alpha_loading, betas_thresholds, irt, mode = "find_intersection"),
        betas_thresholds,
        rtol = rtol,
        atol = atol,
        ctol = ctol
      )$root
    } else if (feature == "half") {
      rootSolve::multiroot(
        gpcm_half_fun(alpha_loading, betas_thresholds, irt, mode = "find_parameters"),
        betas_thresholds,
        rtol = rtol,
        atol = atol,
        ctol = ctol
      )$root
    }
  }
}

convert_pars_list <- function(pars_list, to) {
  conv_std <- function(alphas, betas, mean_ability, sd_ability, to) {
    # if (parameterization = "fa") {
    #   if (center == "thresholds") {
    #
    #   }
    #   if (lv == "center") {
    #     mean_ability_new <- 0.0
    #     sd_ability_new <- sd_ability
    #     betas_new <- betas - mean_ability
    #   }
    # }

    if (to == "bh") {
      list(
        mean_ability = (mean_ability - mean(betas)) * alphas[1],
        sd_ability = sd_ability * alphas[1],
        alphas = alphas / alphas[1],
        betas = (betas - mean(betas)) * alphas[1]
      )
    } else if (to == "stn") {
      list(
        mean_ability = 0.0,
        sd_ability = 1.0,
        alphas = alphas * sd_ability,
        betas = (betas - mean_ability) / sd_ability
      )
    } else if (to == "atbh") {
      mystery_x <- sum(sweep(betas / sd_ability, 1, alphas, FUN = "*")) / (dim(betas)[2] * sum(alphas))

      list(
        mean_ability = mean_ability / sd_ability - mystery_x,
        sd_ability = 1.0,
        alphas = alphas * sd_ability,
        betas = betas / sd_ability - mystery_x,
        alpha_times_betas = sweep(betas / sd_ability - mystery_x, 1, alphas * sd_ability, "*")
      )
    }
  }

  lapply(
    pars_list,
    function(group)
      append(
        group[c("label", "n")],
        do.call(
          conv_std,
          args = append(
            group[c("alphas", "betas", "mean_ability", "sd_ability")],
            list(to = to)
          )
        )
      )
  )
}

convert_pars_list_to_other_model <- function(pars_list, from, feature, irt) {
  lapply(
    pars_list,
    function(group)
      append(
        group[c("label", "n", "mean_ability", "sd_ability", "alphas")],
        list(
          betas = t(
            sapply(
              1:length(group$alphas),
              function(k)
                convert_pars_to_other_model(
                  alpha_loading = group$alphas[k],
                  betas_thresholds = group$betas[k, ],
                  irt = irt,
                  from = from,
                  feature = feature
                )
            )
          )
        )
      )
  )
}
