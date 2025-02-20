# 
# Obiba's Opal - Upload Testing Datasets
#

library(DSOpal)
library(opalr)
library(tibble)

upload_testing_dataset_table <- function(opal, project_name, table_name, local_file_path) {
    if (! opal.project_exists(opal, project_name))
        opal.project_create(opal, project_name, database = "mongodb")
  
    dataset_name <- load(file = local_file_path)
    dataset      <- eval(as.symbol(dataset_name))
    data         <- as_tibble(dataset, rownames = '_row_id_')
  
    opal.table_save(opal, data, project_name, table_name, id.name = "_row_id_", force = TRUE)
}

opal <- opal.login('username','password', url='url', opts = list(ssl_verifyhost=0, ssl_verifypeer=0))

upload_testing_dataset_table(opal, 'test_gamlss_project', 'gamlss1', testthat::test_path('data_files', 'GAMLSS', 'gamlss1.rda'))
upload_testing_dataset_table(opal, 'test_gamlss_project', 'gamlss2', testthat::test_path('data_files', 'GAMLSS', 'gamlss2.rda'))
upload_testing_dataset_table(opal, 'test_gamlss_project', 'gamlss3', testthat::test_path('data_files', 'GAMLSS', 'gamlss3.rda'))
upload_testing_dataset_table(opal, 'test_gamlss_project', 'gamlss_red', testthat::test_path('data_files', 'GAMLSS', 'gamlss_red.rda'))
upload_testing_dataset_table(opal, 'test_gamlss_project', 'gamlss_na', testthat::test_path('data_files', 'GAMLSS', 'gamlss_na.rda'))

opal.logout(opal)
