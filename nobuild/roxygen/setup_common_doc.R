
classes <- list(list(name = "SampleForestNode",
                     type = "forest_node"),
                list(name = "PhylogeneticForestNode",
                     type = "forest_node"),
                list(name = "SampleForest",
                     type = "forest"),
                list(name = "PhylogeneticForest",
                     type = "forest"),
                list(name = "SampleForestNodeTour",
                     type = "node_tour"),
                list(name = "PhylogeneticForestNodeTour",
                     type = "node_tour"))

disabled <- list()

template_dir <- "nobuild/roxygen/txt/"
R_dir <- "R"

header_filename <- file.path(template_dir, "copyright_header.txt")
header <- readLines(header_filename, warn = FALSE)

for (class_list in classes) {
  if (!class_list$type %in% disabled) {
    output_filename <- file.path(R_dir,
                                 paste0(class_list$type, "_doc.R"))

    if (file.exists(output_filename)) {
      file.remove(output_filename)
    }

    writeLines(header, output_filename)
  }
}

for (class_list in classes) {
  if (!class_list$type %in% disabled) {
    input_filename <- file.path(template_dir,
                                paste0(class_list$type, "_doc.txt"))

    template <- readLines(input_filename, warn = FALSE)

    output <- gsub("${CLASS_TYPE}", class_list$name, template, fixed = TRUE)

    output_filename <- file.path(R_dir, paste0(class_list$type, "_doc.R"))

    cat(output, file = output_filename, sep = "\n", append = TRUE)
  }
}
