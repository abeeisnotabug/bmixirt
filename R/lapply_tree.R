make_lapply_tree_expr <- function(iterators, iterator_names, body) {
  if (length(iterators) > 1) {
    sprintf(
      "lapply(%s, function(%s) {%s})",
      iterators[1],
      iterator_names[1],
      make_lapply_tree_expr(iterators[-1], iterator_names[-1], body)
    )
  } else {
    sprintf(
      "lapply(%s, function(%s) %s)",
      iterators[1],
      iterator_names[1],
      body
    )
  }
}

lapply_tree <- function(iterators, iterator_names, my_expr) {
  if (!(is.character(iterators) & is.character(iterator_names))) {
    stop("iterators and iterator_names must be character vectors.")
  } else if (length(iterators) != length(iterator_names)) {
    stop("iterators and iterator_names must have same length.")
  } else {
    eval(
      str2expression(
        make_lapply_tree_expr(iterators, iterator_names, body = paste(deparse(substitute(my_expr, env = environment())), collapse = "\n"))
      ),
      envir = parent.frame(n = 2)
    )
  }
}
