# Task guide (nm_guide): every function and every argument cited in the
# examples of inst/guide/guide_tasks.csv must exist, so that the guide does
# not break silently when the package changes.

guide_tasks <- function() {
  f <- system.file("guide", "guide_tasks.csv", package = "nmpilot")
  utils::read.csv(f, stringsAsFactors = FALSE, encoding = "UTF-8")
}

# all calls of an expression: list of (function name, argument names)
guide_calls <- function(e) {
  out <- list()
  walk <- function(x) {
    if (is.call(x)) {
      fn <- x[[1]]
      if (is.name(fn)) {
        nm <- names(as.list(x))[-1]
        out[[length(out) + 1L]] <<- list(fun = as.character(fn),
                                          args = nm[!is.na(nm) & nzchar(nm)])
      }
      xs <- as.list(x)
      for (k in seq_along(xs))      # empty arguments (x[, 1]) are skipped
        if (!identical(xs[[k]], quote(expr = ))) walk(xs[[k]])
    } else if (is.expression(x)) {
      for (k in seq_along(x)) walk(x[[k]])
    }
  }
  walk(e)
  out
}

ns <- asNamespace("nmpilot")
exports <- getNamespaceExports("nmpilot")

# formals accepted for a call: the function itself, or the union of the
# nmpilot methods for the generics plot, summary and print
guide_formals <- function(fun) {
  if (fun %in% c("plot", "summary", "print")) {
    m <- ls(ns, all.names = TRUE, pattern = paste0("^", fun, "\\.nm_"))
    return(unique(unlist(lapply(m, function(z) names(formals(get(z, envir = ns)))))))
  }
  names(formals(get(fun, envir = ns)))
}

test_that("the guide file is well formed", {
  g <- guide_tasks()
  expect_true(all(c("task_id", "stage_id", "stage", "task", "functions",
                    "description", "code", "status") %in% names(g)))
  expect_false(any(duplicated(g$task_id)))
  expect_true(all(g$status %in% c("available", "planned")))
  av <- g$status == "available"
  expect_true(all(nzchar(g$code[av])))
  expect_true(all(nzchar(g$functions[av])))
})

test_that("functions named by the available tasks exist in nmpilot", {
  g <- guide_tasks()
  g <- g[g$status == "available", ]
  f <- unique(trimws(unlist(strsplit(g$functions, ",", fixed = TRUE))))
  f <- f[nzchar(f)]
  missing <- f[!vapply(f, function(z) z %in% exports ||
                         exists(z, envir = ns, inherits = FALSE), logical(1))]
  expect_equal(missing, character(0))
})

test_that("every example parses, and its nmpilot calls use existing arguments", {
  g <- guide_tasks()
  g <- g[g$status == "available", ]
  problems <- character(0)
  for (i in seq_len(nrow(g))) {
    e <- tryCatch(parse(text = g$code[i], keep.source = FALSE),
                  error = function(err) err)
    if (inherits(e, "error")) {
      problems <- c(problems, sprintf("task %s: code does not parse", g$task_id[i]))
      next
    }
    for (cl in guide_calls(e)) {
      is_nm <- startsWith(cl$fun, "nm_") || cl$fun %in% c("plot", "summary", "print")
      if (!is_nm) next
      if (startsWith(cl$fun, "nm_") && !(cl$fun %in% exports)) {
        problems <- c(problems, sprintf("task %s: %s() is not exported",
                                        g$task_id[i], cl$fun))
        next
      }
      ok <- guide_formals(cl$fun)
      # an argument caught by `...` is a former (deprecated) name or a typo
      bad <- setdiff(cl$args, setdiff(ok, "..."))
      if (length(bad))
        problems <- c(problems, sprintf("task %s: %s(%s) unknown argument(s)",
                                        g$task_id[i], cl$fun,
                                        paste(bad, collapse = ", ")))
    }
  }
  expect_equal(problems, character(0))
})

test_that("nm_guide writes the page with every task", {
  f <- tempfile(fileext = ".html")
  nm_guide(output_file = f, open = FALSE)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  g <- guide_tasks()
  expect_equal(lengths(regmatches(h, gregexpr("<article class=\"detail\"", h))), nrow(g))
})
