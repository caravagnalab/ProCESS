devtools::load_all()

generate_sample_forest_example <- function() {
  # set the seed of the random number generator
  set.seed(0)

  sim <- TissueSimulation("Two Populations", epigenetic_states = c("E1", "E2"))

  sim$death_activation_level <- 20

  # add mutant "A" and set its species rates
  sim$add_mutant("A", list(E1 = list(duplication = 0.1, death = 0.1,
                                     E2 = 0.01),
                           E2 = list(duplication = 0.08, death = 0.01,
                                     E1 = 0.01)))

  sim$place_cell("A[E1]", 500, 500)
  sim$run_up_to_size("A[E1]", 1000)

  bbox_width <- 10

  sim$sample_cells("S_1_1",
                   bottom_left = c(480, 480),
                   top_right = c(480 + bbox_width, 480 + bbox_width))

  sim$sample_cells("S_1_2",
                   bottom_left = c(500, 500),
                   top_right = c(500 + bbox_width, 500 + bbox_width))

  sim$run_up_to_time(sim$get_clock() + 15)

  cell <- sim$choose_border_cell_in("A")

  sim$add_mutant("B", list(E1 = list(duplication = 0.08, death = 0.05,
                                     E2 = 0.05),
                           E2 = list(duplication = 0.3, death = 0.05,
                                     E1 = 0.1)))

  sim$mutate_progeny(cell, "B")

  # let it grow more
  sim$run_up_to_size("B[E1]", 7000)

  n_w <- n_h <- 25
  ncells <- 0.9 * n_w * n_h

  bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
  sim$sample_cells("S_2_1", bbox$lower_corner, bbox$upper_corner)

  bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
  sim$sample_cells("S_2_2", bbox$lower_corner, bbox$upper_corner)

  sim$get_sample_forest()
}

generate_phylo_forest_example <- function(sample_forest) {
  # set the seed of the random number generator
  set.seed(0)

  # building a mutation engine by using the "demo" setup
  m_engine <- MutationEngine(setup_code = "demo")

  m_engine$add_mutant(mutant_name = "A",
                      passenger_rates = list("E1" = c(SNV = 1e-9, indel = 1e-9),
                                             "E2" = c(SNV = 3e-9, CNA = 1e-11)),
                      drivers = list(SNV("22", 46510210, "C", allele = 1),
                                     Mutation("22", 16085675, "GCCTCCCGA", "G"),
                                     CNA(type = "A", chr = "22",
                                         from = 10303470, len = 200000),
                                     CNA("D", "22", 5010000, 200000,
                                         allele = 1),
                                     WGD))

  m_engine$add_mutant("B", list("E1" = c(SNV = 8e-9), "E2" = c(SNV = 2e-8)),
                      list(list("DGCR8 P26L", allele = 1),  # the first SNV
                           SNV("22", 12028576, "G")))       # the second SNV

  m_engine$add_exposure(coefficients = c(SBS13 = 0.2, SBS1 = 0.8))
  m_engine$add_exposure(c(ID2 = 0.6, ID13 = 0.2, ID21 = 0.2))

  m_engine$add_exposure(time = 100, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))

  m_engine$add_exposure(time = 120, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5,
                                      ID1 = 0.8, ID9 = 0.2))

  subjects <- m_engine$get_germline_subjects()
  m_engine$set_germline_subject(subjects[2, "sample"])

  m_engine$place_mutations(sample_forest, 1000, 500)
}

sample_forest <- generate_sample_forest_example()
phylo_forest <- generate_phylo_forest_example(sample_forest)

# create the data dir
data_dir <- file.path("inst", "extdata")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

sample_forest$save(file.path(data_dir, "sample_forest_example.sff"))
phylo_forest$save(file.path(data_dir, "phylo_forest_example.pff"))

# move reference and reference index from "demo" to inst/extdata
reference_path <- file.path(data_dir, "example_ref.fasta")
file.copy(phylo_forest$get_reference_path(), reference_path)
reference_chi_path <- paste0(reference_path, ".chi")
file.copy(paste0(phylo_forest$get_reference_path(), ".chi"),
          paste0(reference_path, ".chi"))

library(R.utils)

gzip(reference_path)
