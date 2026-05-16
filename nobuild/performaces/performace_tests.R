tests <- c("mutation_engine-time.R", "mutants-time.R", "num_of_cells-time.R",
           "sample_forest-time.R", "phylo_forest-time.R", "sequencing-time.R")

if (Sys.getenv("RSTUDIO") == "1") {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  script_dir <- getwd()
}

image_dir <- file.path(script_dir, "images")
working_dir <- file.path(script_dir, "working_dir")
result_dir <- file.path(script_dir, "results")
test_dir <- file.path(script_dir, "tests")

for (directory_name in c(image_dir, working_dir, result_dir)) {
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
  }
}

for (test in tests) {
  vars <- ls()
  rm(list=vars[!(vars %in% c("test", "tests", 
                             "image_dir", "working_dir", "result_dir",
                             "test_dir"))])
  setwd(test_dir)
  
  message(paste("Running", test, "..."))
  source(test)
}
