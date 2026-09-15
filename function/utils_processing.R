# utils_processing.R for data processing functions

# Function to update status and progress bar
update_status <- function(session, messages, new_message, progress_value) {
  print(new_message)
  status(new_message)
  updateProgressBar(session, id = "progress", value = progress_value)
  messages(c(messages(), new_message))
  Sys.sleep(1)  # Simulate delay for UI responsiveness
}

# Function to handle each step with error handling
run_step <- function(step_function, session, progress_value, ...) {
  tryCatch({
    step_function(...)
    updateProgressBar(session, id = "progress", value = progress_value)
  }, error = function(e) {
    error_message <- paste("Error:", e$message)
    messages(c(messages(), error_message))
  })
}
