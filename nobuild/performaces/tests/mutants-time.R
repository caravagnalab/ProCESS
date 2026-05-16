devtools::load_all()
library(peakRAM)
library(dplyr)
library(numbers)
library(ggplot2)

source("global.R")

setup_tumour <- function(border_growth_model)
{
  width <- 7000
  height <- 7000
  
  sim <- TissueSimulation("Test", width, height)
  
  sim$border_growth_model <- border_growth_model
  
  return(sim)
}

setup_multi_clone_tumour <- function(num_of_clones,
                                     border_growth_model, grow_rate,
                                     num_of_slices)
{
  sim <- setup_tumour(border_growth_model)

  tissue_size <- sim$get_tissue_size()

  alive_formula <- NULL
  for (s in 1:num_of_slices) {
    i <- s %% num_of_clones;
    mutant_name <- paste("Clone", i)

    x <- as.integer((cos(s*pi/num_of_slices)/2+1)*tissue_size[1]/2)
    y <- as.integer((sin(s*pi/num_of_slices)/2+1)*tissue_size[1]/2)

    if (nrow(sim$get_added_cells() %>% filter(mutant==mutant_name))==0) {
      sim$add_mutant(name = mutant_name,
                     growth_rates = grow_rate, death_rates = 0.01)
      
      if (is.null(alive_formula)) {
        alive_formula <- sim$var(mutant_name)
        overall_formula <- 2*sim$var(paste0(mutant_name,".duplications")) + 1
      } else {
        alive_formula <- alive_formula + sim$var(mutant_name)
        overall_formula <- overall_formula +
          2*sim$var(paste0(mutant_name,".duplications")) + 1
      }
    }

    sim$place_cell(mutant_name, x, y)
  }

  return(list("sim"=sim, "alive_formula"=alive_formula,
              "overall_formula"=overall_formula))
}

grow_multi_clone_clones <- function(sim, formula, size_goal)
{
  sim$run_until(formula >= size_goal)
  
  return(sim)
}

complete_grow_multi_clone_clones <- function(formula, size_goal,
                                             initial_size, num_of_clones,
                                             border_growth_model,
                                             init_growth_rate, growth_rate)
{
  result <- setup_multi_clone_tumour(num_of_clones,
                                     border_growth_model, growth_rate)
  
  return(grow_multi_clone_clones(result$sim, result[[formula]], size_goal))
}

add_new_data <- function(sim, time, growth_rate_v,
                         border_growth_v, repetition_id_v,
                         num_of_clones_v, output_file, total_result=NULL)
{
  counts <- sim$get_counts()
  
  num_of_alive_clones <- sim$get_counts() %>% filter(counts >0) %>% nrow()
  
  result <- data.frame(elapsed_time=time, growth_rate=growth_rate_v,
                       num_of_alive_cells=sum(counts$counts),
                       num_of_overall_cells=sum(counts$overall),
                       border_growth_model=border_growth_v,
                       repetition_id=repetition_id_v,
                       num_of_clones=num_of_clones_v,
                       alive_clones=num_of_alive_clones,
                       simulated_time=sim$get_clock())
  
  if (is.null(total_result)) {
    total_result <- result
  } else {
    total_result <- rbind(total_result, result)
  }
  
  saveRDS(total_result, output_file)
  
  return(total_result)
}

final_size <- 1000000
num_of_repetitions <- 5
border_growth_v <- TRUE
growth_rate_v <- 0.4

set.seed(0)

num_of_mutant_list <- c(10, 20, 30, 40, 48, 60, 80, 120, 240)
repetition_list <- 1:num_of_repetitions

test_df <- expand.grid(num_of_mutant_list, repetition_list)

if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

output_file <- file.path(result_dir, "mutants-time.rds")

section <- "multi-mutant"
img_ext <- "png"

total_result <- NULL

if (file.exists(output_file)) {
  total_result <- readRDS(output_file)
}

for (test_idx in sample(1:nrow(test_df))) {
  num_of_mutants <- test_df[test_idx, 1]
  
  j <- test_df[test_idx, 2]
  
  to_measure <- is.null(total_result)
  
  if (!to_measure) {
    num_of_measures <- nrow(total_result %>%
                              filter(growth_rate == growth_rate_v &
                                       border_growth_model == border_growth_v &
                                       num_of_clones == num_of_mutants &
                                       repetition_id == j))
    
    to_measure <- (num_of_measures == 0)
  }
  
  if (to_measure) {
    set.seed(j)
    cat("Setup:\n")
    cat(paste("  num_of_mutants:", num_of_mutants, "\n"))
    cat(paste("  num_of_mutants:", num_of_mutants, "\n"))
    cat(paste("  repetition_id:", j, "\n"))
      setup <- setup_multi_clone_tumour(num_of_mutants, border_growth_v,
                                        growth_rate_v,
                                        mLCM(num_of_mutant_list))

      formula <- setup[["alive_formula"]]
      
      result <- peakRAM(sim <- grow_multi_clone_clones(setup$sim, formula,
                                                       final_size))

      total_result <- add_new_data(sim, result$Elapsed_Time_sec,
                                   growth_rate_v,
                                   border_growth_v, j,
                                   num_of_mutants, output_file, total_result)
  }
}

plot <- ggplot(data = total_result,
               mapping = aes(x = alive_clones, y = elapsed_time/60)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 7) +
  stat_summary(fun = "mean", geom = "point", 
               size = 3, color = "red") +
  labs(
    x = "Mutants",
    y = "Computation time (minutes)"
  ) +
  geom_smooth(method = "lm", formula = y ~ poly(x,1), color = "black") +
  theme(legend.position = "bottom") +
  theme_minimal()

save_image(section, "mutants-vs-time", extension=img_ext,
           directory=image_dir)

results <- compare_complexity(total_result,
                              "alive_clones", "elapsed_time",
                              "Alive clones-vs-Computation time")

plot2 <- ggplot(total_result, aes(x = alive_clones,
                                  y = elapsed_time/60)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2) +
  stat_summary(fun = "mean",geom = "point", size = 3, color = "red") +
  labs(
    x = "Mutants",
    y = "Computation time (minutes)"
  ) +
  theme_minimal()

results <- compare_complexity(total_result,
                              "num_of_clones", "elapsed_time",
                              "Num. of clones-vs-Computation time")

cat("\n")
