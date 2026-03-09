suppressPackageStartupMessages({
  library(bench)
})

benchmark_wrapper <- function(expr, iterations = 3, envir = parent.frame()) {
  dur <- numeric(iterations)
  cpu_percent <- numeric(iterations)

  expr <- base::substitute(expr)

  for (i in seq_len(iterations)) {
    tm <- system.time(eval(expr, envir = envir))
    dur[i] <- tm[["elapsed"]]
    cpu_percent[i] <- (tm[["user.self"]] + tm[["sys.self"]]) / tm[["elapsed"]]
  }

  res <- bench::mark(
    eval(expr, envir = envir),
    iterations = 1,
    check = FALSE,
    memory = TRUE,
    time_unit = "s"
  )

  data.frame(
    duration = median(dur),
    cpu_usage = median(cpu_percent),
    memory = as.numeric(res$mem_alloc)
  )
}
