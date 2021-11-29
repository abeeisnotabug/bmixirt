make_lapply_tree_expr <- function(iterators, body) {
  if (length(iterators) > 1) {
    sprintf(
      "lapply(%s, function(%s) {%s})",
      iterators[1],
      names(iterators[1]),
      make_lapply_tree_expr(iterators[-1], body)
    )
  } else {
    sprintf(
      "lapply(%s, function(%s) %s)",
      iterators[1],
      names(iterators[1]),
      body
    )
  }
}

lapply_tree <- function(iterators, my_expr) {
  if (!is.character(iterators) || is.null(names(iterators))) {
    stop("iterators must be a named character vector.")
  } else {
    eval(
      str2expression(
        make_lapply_tree_expr(iterators, body = paste(deparse(substitute(my_expr, env = environment())), collapse = "\n"))
      ),
      envir = parent.frame(n = 2)
    )
  }
}
