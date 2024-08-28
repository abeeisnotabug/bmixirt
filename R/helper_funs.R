give_pois_mccirt <- function(log = TRUE, cent = TRUE, marginal = FALSE)
  c("uncond_classprob", "intercept", paste0("threshold", if (cent) "_c"),
    paste0(if (log) "log_", "loading_", c("I", "R", "T"), "_c"),
    paste0(if (log) "log_", "scaling_", c("I", "R", "T")),
    "ability_T", if (!marginal) c("ability_R", "ability_I"))

give_pois_mirt <- function(log = TRUE, cent = TRUE, marginal = FALSE)
  c("uncond_classprob", "intercept", paste0("threshold", if (cent) "_c"),
    paste0(if (log) "log_", "loading_", "_c"),
    paste0(if (log) "log_", "scaling_"),
    if (!marginal) "ability")

make_run_model_code <- function(model_name, data_name, stan_file, interface = "cmdstanr", C, n_chains = 4, n_iter = 10000, n_warmup = 4000, refresh = 100, adapt_delta = .8, max_treedepth = 10, init = "0", seed = 1909, integrate = 1, nohup = FALSE, disengage = TRUE)
  sprintf(
    "%sRscript --vanilla run_mixmodel.R %s %s %s %s %i %i %i %i %i %.3f %i %s %i %i > fits/logfiles/%s/%s.out 2> fits/logfiles/%s/%s.err%s",
    if (nohup) "nohup " else "",
    model_name, data_name, stan_file, interface, C, n_chains, n_iter, n_warmup, refresh, adapt_delta, max_treedepth, init, seed, as.integer(integrate),
    interface, model_name,
    interface, model_name,
    if (disengage) " &" else ""
  )

run_model <- function(model_name, data_name, stan_file, interface, C, n_chains, n_iter, n_warmup, refresh, adapt_delta, max_treedepth, init, seed, integrate)
  system(make_run_model_code(model_name, data_name, stan_file, interface, C, n_chains, n_iter, n_warmup, refresh, adapt_delta, max_treedepth, init, seed, integrate, nohup = TRUE, disengage = FALSE), intern = FALSE, wait = FALSE)

extract_itempars <- function(draws) {
  intercept_logscalingR <- as.data.table(posterior::subset_draws(draws, variable = c("intercept", "log_scaling_R")))
  intercept_logscalingR[, c(".chain", ".iteration") := NULL]
  intercept_logscalingR <- melt(intercept_logscalingR, measure = patterns("^intercept", "^log_scaling_R"),
                                variable.name = "c", value.name = c("intercept", "log_scaling_R"))

  logscalingIT <- as.data.table(posterior::subset_draws(draws, variable = c("log_scaling_I", "log_scaling_T")))
  logscalingIT[, c(".chain", ".iteration") := NULL]

  thresholdsc <- as.data.table(posterior::subset_draws(draws, variable = "threshold_c"))
  thresholdsc[, c(".chain", ".iteration") := NULL]
  setnames(thresholdsc, names(thresholdsc), str_replace_all(names(thresholdsc), c("\\[" = ",", "\\]" = "")))
  thresholdsc <- melt(thresholdsc, measure = patterns("^threshold_c"), value.name = "threshold_c")
  thresholdsc[, c("variable", "c", "item", "category") := tstrsplit(variable, ",", fixed = TRUE, type.convert = FALSE)]
  thresholdsc[, c("c", "item") := .(factor(c), factor(item, levels = seq_along(unique(item)), ordered = TRUE))]
  thresholdsc <- dcast(thresholdsc, .draw + c + item ~ variable + category, value.var = "threshold_c")

  logloadingsc <- as.data.table(posterior::subset_draws(draws, variable = c("log_loading_I_c", "log_loading_R_c", "log_loading_T_c")))
  logloadingsc[, c(".chain", ".iteration") := NULL]
  setnames(logloadingsc, names(logloadingsc), str_replace_all(names(logloadingsc), c("\\[" = ",", "\\]" = "")))
  logloadingsc <- melt(logloadingsc, measure = patterns("^log_loading_I_c", "^log_loading_R_c", "^log_loading_T_c"),
                       variable.name = "item", value.name = c("log_loading_I_c", "log_loading_R_c", "log_loading_T_c"))
  logloadingsc[, item := factor(item, ordered = TRUE)]

  ucp <- as.data.table(posterior::subset_draws(draws, variable = "uncond_classprob"))
  ucp[, c(".chain", ".iteration") := NULL]
  ucp <- melt(ucp, measure = patterns("^uncond_classprob"), variable.name = "c", value.name = "uncond_classprob")
  ucp[, c := factor(c, labels = seq_along(unique(c)))]

  parDT <- merge(intercept_logscalingR, logscalingIT)
  parDT <- merge(parDT, logscalingIT)
  parDT <- merge(parDT, thresholdsc, by = c(".draw", "c"))
  parDT <- merge(parDT, logloadingsc, by = c(".draw", "item"))
  parDT <- merge(parDT, ucp, by = c(".draw", "c"))

  scaling_names <- str_sort(str_subset(names(parDT), "log_"))
  threshold_names <- str_sort(str_subset(names(parDT), "threshold"), numeric = TRUE)
  setcolorder(parDT, c(".draw", "c", "item", rev(scaling_names), "intercept", threshold_names, "uncond_classprob"))
  # parDT[, c := dplyr::recode(c, !!!setNames(order(.SD[, var(c(as.matrix(.SD))), .SDcols = threshold_names, by = c]$V1), seq_along(unique(c)))), by = .draw]
  setkey(parDT, .draw, c, item)
  parDT[]
}

plot_iccs <- function(itempars, xmin, xmax, xstep = .01, filtermodel = c("grm", "gpcm")) {
  itempars <- itempars[, .(x = seq(xmin, xmax, xstep)), by = eval(colnames(itempars))]
  setkey(itempars, model, c, item)

  threshold_names <- str_subset(colnames(itempars), "^threshold")

  itempars[, cloading := cloading * x]
  itempars[, (threshold_names) := lapply(.SD, function(th) cloading - th), .SDcols = threshold_names]
  itempars[, (str_subset(colnames(itempars), "loading")) := NULL]

  threshold_0_name <- str_replace(threshold_names[1], "1", "0")
  itempars["gpcm", (threshold_0_name) := 0]
  threshold_names <- str_sort(str_subset(colnames(itempars), "^threshold"), numeric = TRUE)
  setcolorder(itempars, c("model", "c", "item", threshold_names))

  itempars["gpcm", (threshold_names) := purrr::accumulate(.SD, `+`), .SDcols = threshold_names]
  itempars["gpcm", (threshold_names) := exp(.SD), .SDcols = threshold_names]
  itempars["gpcm", (threshold_names) := .SD / rowSums(.SD), .SDcols = threshold_names]

  itempars["grm", (threshold_names[-1]) := 1 / (1 + exp(-.SD)), .SDcols = threshold_names[-1]]
  itempars["grm", (threshold_0_name) := 1]
  itempars["grm", head(threshold_names, -1) := subset(.SD, select = head(threshold_names, -1)) - subset(.SD, select = tail(threshold_names, -1)), .SDcols = threshold_names]

  itempars <- melt(itempars, measure = patterns("^threshold"), variable.name = "category")

  ggplot(itempars[model %in% filtermodel], aes(x = x, y = value, color = model, group = interaction(model, category))) +
    geom_line() +
    geom_function(fun = dnorm, color = "black", inherit.aes = FALSE) +
    scale_color_manual(values = c(rgb(1.0, 0.5, 0.0), rgb(0.0, 0.5, 0.5))) +
    facet_grid(rows = vars(c), cols = vars(item))
}

grm_diff_fun <- function(intersections_to_find, parms) {
  Jp1 <- length(parms)
  J <- Jp1 - 1

  pred <- parms[1] * outer(intersections_to_find, parms[2:Jp1], `-`)

  logits <- cbind(1, 1 / (1 + exp(-pred)))
  category_probs <- cbind(logits[, 1:J] - logits[, 2:Jp1], logits[, Jp1])
  differences <- category_probs[, 1:J] - category_probs[, 2:Jp1]

  diag(differences)
}

extract_reorder_draws <- function(draws, metadata, ordercols = "var_threshold_c", reorder = TRUE) {
  cat(sprintf("## Reordering parameters based on %s: ", paste(ordercols, collapse = ", ")), sep = "")

  id.vars <- c(".chain", ".iteration", ".draw")
  vars_in_draws <- unique(str_extract(posterior::variables(draws), "[^\\[]+"))

  corder <- as.data.table(posterior::subset_draws(draws, variable = ordercols))
  corder <- melt(corder, measure.vars = patterns(ordercols), variable.name = "c", value.name = ordercols)
  corder[, c := as.character(factor(c, labels = seq_along(levels(c))))]

  if (length(ordercols) == 1) {
    corder[, c_re := as.character(order(order(get(ordercols)))), by = id.vars]
  } else if (length(ordercols) == 2) {
    corder[, c_re := as.character(order(c(which.min(get(ordercols[1])), setdiff(order(get(ordercols[2])), which.min(get(ordercols[1])))))), by = id.vars]
  } else {
    stop("Sorry, more than two ordercols are not supported because the re-ordering is hard-coded lol.")
  }

  mixpars_1 <- names(which(sapply(metadata$stan_variable_sizes, \(sizevec) sizevec[1] == metadata$stan_variable_sizes$uncond_classprob)))
  mixpars_2 <- names(which(sapply(metadata$stan_variable_sizes, \(sizevec) sizevec[2] == metadata$stan_variable_sizes$uncond_classprob)))
  nomixpars <- setdiff(names(metadata$stan_variable_sizes), c(mixpars_1, mixpars_2))

  if (!setequal(c(mixpars_1, mixpars_2, nomixpars), metadata$stan_variables))
    stop("!!!!!!!!!! The parameters found do not match.")

  mixpars_1 <- intersect(mixpars_1, vars_in_draws)
  mixpars_2 <- intersect(mixpars_2, vars_in_draws)
  nomixpars <- intersect(nomixpars, vars_in_draws)

  if (!length(c(mixpars_1, mixpars_2)))
    stop("!!!!!!!!!! There are no mixture parameters in the draws.")

  if (!isTRUE(all.equal(corder[, uniqueN(c_re), by = .(.chain, c)][, unique(V1)], 1)))
    warning("!!!!! There seems to be more than one ordering per group.")

  if (reorder) {
    combdraws <- as.data.table(posterior::subset_draws(draws, variable = nomixpars))

    for (parname in c(mixpars_1, mixpars_2)) {
      cat(paste0(parname, if (parname != tail(c(mixpars_1, mixpars_2), 1)) ", "))

      this_DT <- as.data.table(posterior::subset_draws(draws, variable = parname))
      this_DT <- melt(this_DT, id.vars = id.vars)
      this_DT[, idx := variable]
      this_DT[, idx := as.character(factor(idx, labels = str_sub(str_split_fixed(levels(idx), "\\[", 2)[, 2], 1, -2)))]
      this_DT[, variable := as.character(factor(variable, labels = str_split_fixed(levels(variable), "\\[", 2)[, 1]))]
      this_DT[, c := as.character(idx)]
      if (parname %in% mixpars_1) {
        Cl <- str_length(metadata$stan_variable_sizes$uncond_classprob)
        this_DT[, c := str_sub(idx, 1, Cl)]
        this_DT[, idx := str_sub(idx, Cl + 1, -1)]
        this_DT[corder, c := c_re, on = c(".chain", ".iteration", ".draw", "c")]
        this_DT[, variable := paste0(variable, "[", c, idx, "]")]
      } else {
        this_DT[, c := str_sub(idx, -Cl, -1)]
        this_DT[, idx := str_sub(idx, 1, -(Cl + 1))]
        this_DT[corder, c := c_re, on = c(".chain", ".iteration", ".draw", "c")]
        this_DT[, variable := paste0(variable, "[", idx, c, "]")]
      }
      this_DT <- dcast(this_DT, .chain + .iteration + .draw ~ variable, value.var = "value")

      combdraws[this_DT, (setdiff(names(this_DT), id.vars)) := mget(setdiff(names(this_DT), id.vars)), on = id.vars]
    }
    cat("\n## Done.", sep = "\n")

    as_draws_df(combdraws)
  }
}
