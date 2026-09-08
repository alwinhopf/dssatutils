# Shared bounded recovery for transient service/network failures. Exhausting
# the budget stops the provider stream; it must not mark the remaining points
# as missing geographic coverage.
.provider_transient <- function(error) {
  grepl("resolve host|resolv.*host|name resolution|nodename nor servname|connection|connect to|timed? ?out|timeout|429|50[234]|temporarily unavailable",
        conditionMessage(error), ignore.case = TRUE)
}
.provider_retry <- function(operation, attempts = 3L, delay_seconds = 5, sleep = Sys.sleep) {
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(list(value = operation()), error = function(e) e)
    if (!inherits(result, "error")) return(result$value)
    if (!.provider_transient(result)) stop(result)
    if (attempt == attempts) {
      stop(structure(list(message = paste("Provider connectivity unavailable after bounded retries:",
                                        conditionMessage(result)), call = NULL),
                     class = c("dssat_connectivity_error", "error", "condition")))
    }
    wait <- min(delay_seconds * 2^(attempt - 1L), 30)
    message(sprintf("Provider connection failed; retry %d/%d in %.0f seconds", attempt + 1L, attempts, wait))
    sleep(wait)
  }
}
