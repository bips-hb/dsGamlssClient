setupGAMLSSTest <- function(packages=c(), env=parent.frame())
{
  datasets <- c("gamlss1", "gamlss2", "gamlss3")
  logindata <- "logindata.dslite.gamlss"
  dslite.server <- "dslite.server"

  # check server-side package is installed
  for (package in packages) {
    if (!base::requireNamespace(package, quietly = TRUE)) {
      stop(package, " package is required in the local R installation for the execution of the tests.",
           call. = FALSE)
    }
  }
  
  # load simulated test datasets and corresponding login definition
  load(testthat::test_path("data_files", "GAMLSS", paste(datasets[1], ".rda", sep="")), envir=env)
  load(testthat::test_path("data_files", "GAMLSS", paste(datasets[2], ".rda", sep="")), envir=env)
  load(testthat::test_path("data_files", "GAMLSS", paste(datasets[3], ".rda", sep="")), envir=env)
  load(testthat::test_path("data_files", "GAMLSS", paste(logindata, ".rda", sep="")), envir=env)
  
  # new DSLiteServer, hosting the simulated test datasets
  tables <- list()
  for (dataset in datasets) {
    tables[[dataset]] <- get(dataset, envir = env)
  }
  
  # get server symbol from login data if not provided by argument
  dslite.server.sym <- dslite.server
  if (is.null(dslite.server)) {
    ld <- get(logindata, envir = env)
    dslite.server.sym <- as.character(ld[!is.null(ld$url),][1,]$url)
  }
  
  assign(dslite.server, newDSLiteServer(tables=tables, config = DSLite::defaultDSConfiguration(include=packages)), envir = env)
  # specify in which environment the dslite server lives
  options(datashield.env=env)
  
  # return logindata
  get(logindata, envir=env)
  
}

init.studies.dataset.gamlss <- function(variables)
{
  if (ds.test_env$secure_login_details)
  {
    #reading data from local files
    ds.test_env$local.values.1 <- read.csv(testthat::test_path("data_files", "GAMLSS", "gamlss1.csv"), header = TRUE)
    ds.test_env$local.values.2 <- read.csv(testthat::test_path("data_files", "GAMLSS", "gamlss2.csv"), header = TRUE)
    ds.test_env$local.values.3 <- read.csv(testthat::test_path("data_files", "GAMLSS", "gamlss3.csv"), header = TRUE)
    ds.test_env$local.values   <- rbind(ds.test_env$local.values.1,ds.test_env$local.values.2,ds.test_env$local.values.3)
    if (ds.test_env$driver == "OpalDriver")
    {
      builder <- DSI::newDSLoginBuilder(.silent = TRUE)
      builder$append(server = "server1", url = ds.test_env$ip_address_1, user = ds.test_env$user_1, password = ds.test_env$password_1, profile = ds.test_env$profile_1, table = "test_gamlss_project.gamlss1", options=ds.test_env$options_1)
      builder$append(server = "server2", url = ds.test_env$ip_address_2, user = ds.test_env$user_2, password = ds.test_env$password_2, profile = ds.test_env$profile_2, table = "test_gamlss_project.gamlss2", options=ds.test_env$options_2)
      builder$append(server = "server3", url = ds.test_env$ip_address_3, user = ds.test_env$user_3, password = ds.test_env$password_3, profile = ds.test_env$profile_3, table = "test_gamlss_project.gamlss3", options=ds.test_env$options_3)
      ds.test_env$login.data <- builder$build()
    }
    else 
    {
      ds.test_env$login.data <- setupGAMLSSTest(packages=c("dsBase", "dsGamlss", "gamlss", "gamlss.dist"),
                                                env=ds.test_env)
    }
    ds.test_env$stats.var <- variables
    
  }
}

log.in.data.server <- function()
{
  ds.test_env$connections <- datashield.login(logins=ds.test_env$login.data, assign=TRUE, variables=ds.test_env$stats.var, opts=getOption("datashield.opts"))
}

log.out.data.server <- function()
{
  if (!is.null(ds.test_env) && !is.null(ds.test_env$connections))
  {
    datashield.logout(ds.test_env$connections)
  }
  rm(list = ls())
  gc()
}

connect.studies.dataset.gamlss <- function(variables)
{
  log.out.data.server()
  source(testthat::test_path("connection_to_datasets", "login_details.R"))
  init.studies.dataset.gamlss(variables)
  log.in.data.server()
}

disconnect.studies.dataset.gamlss <- function()
{
  log.out.data.server()
}
