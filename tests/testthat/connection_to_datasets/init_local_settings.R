# this file stores some settings for the continuous integration and local testing.

init.ip.address <- function() {
  file.name <- init.local.settings()
  if (file.exists(file.name)) {
    content <- read.csv(file.name, header = FALSE)
    ip.address <- as.character(content[[1]][1])
  } else {
    ip.address <- "localhost"
  }
  return(ip.address)
}



init.local.settings <- function() {
  path <- testthat::test_path("connection_to_datasets", "local_settings.csv")
  return(path)
}
