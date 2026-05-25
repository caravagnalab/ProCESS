devtools::load_all()
library(peakRAM)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(MuMIn)
library(forcats)

source("global.R")

add_new_data <- function(sim, time, growth_rate_v,
                         border_growth_model_v, repetition_id,
                         output_file, total_result=NULL)
{
  counts <- sim$get_counts()
  
  result <- data.frame(elapsed_time=time, growth_rate=growth_rate_v,
                       num_of_alive_cells=counts$counts,
                       num_of_overall_cells=counts$overall,
                       border_growth_model=border_growth_model_v,
                       repetition_id=repetition_id,
                       simulated_time=sim$get_clock())
  
  if (is.null(total_result)) {
    total_result <- result
  } else {
    total_result <- rbind(total_result, result)
  }
  
  saveRDS(total_result, output_file)
  
  return(total_result)
}

if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

output_file <- file.path(result_dir, "num_of_cells-time.rds")

setwd(working_dir)

with_overall <- FALSE
max_alive_cells <- 2000000
max_overall_cells <- 43000000
num_of_samples <- 5
num_of_repetitions <- 5
growth_rates <- c(0.1, 0.2, 0.4)

AICc_delta_threshold <- 1

section <- "tissue"
img_ext <- "png"

total_result <- NULL

if (file.exists(output_file)) {
  total_result <- readRDS(output_file)
}

delta_size <- max_alive_cells / num_of_samples

for (border_growth_v in c(TRUE,FALSE)) {
  for (growth_rate_v in growth_rates) {
    to_measure <- is.null(total_result)

    if (!to_measure) {
      num_of_measures <- nrow(total_result %>%
                                filter(growth_rate == growth_rate_v &
                                         border_growth_model == border_growth_v))
      
      to_measure <- (num_of_measures < num_of_samples*num_of_repetitions)
      
      if (to_measure && num_of_measures>0) {
        total_result <- total_result %>%
          filter(growth_rate != growth_rate_v |
                   border_growth_model != border_growth_v)
        
        saveRDS(total_result, output_file)
      }
    }
    
    if (to_measure) {
      set.seed(0)
      for (repetition_id_v in 1:num_of_repetitions) {
        cat("Setup:\n")
        cat(paste("  border_growth_model:", border_growth_v, "\n"))
        cat(paste("  growth_rate:", growth_rate_v, "\n"))
        cat(paste("  repetition_id:", repetition_id_v, "\n"))
        
        sim <- setup_tumour(border_growth_v, growth_rate_v)
        time <- 0.0
        for (i in 1:num_of_samples) {
          num_of_cells <- as.integer(i*delta_size)
          
          result <- peakRAM(sim$run_up_to_size("Clone 1", num_of_cells))
          
          time <- time + result$Elapsed_Time_sec
          total_result <- add_new_data(sim, time, growth_rate_v,
                                       border_growth_v, repetition_id_v,
                                       output_file, total_result)
        }
  
        if (with_overall) {
          overall_formula <- 2*sim$var("Clone 1.duplications")+1
          result <- peakRAM(sim$run_until(overall_formula>=max_overall_cells))
          time <- time + result$Elapsed_Time_sec
    
          total_result <- add_new_data(sim, time, growth_rate_v,
                                       border_growth_v, repetition_id_v,
                                       output_file, total_result)
        }
      }
    }
  }
}

plot1_data <- total_result %>% filter(growth_rate %in% growth_rates) %>%
  mutate(color_var = interaction(border_growth_model, growth_rate))

plot1_data$color_var <- factor(plot1_data$color_var,
                               levels = c("TRUE.0.1", "TRUE.0.2", "TRUE.0.4",
                                          "FALSE.0.1", "FALSE.0.2",
                                          "FALSE.0.4"))

plot1 <- ggplot(plot1_data,
                  aes(x = num_of_alive_cells, y = elapsed_time/60,
                      color = color_var)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 3e+5) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  labs(
    x = "Tumour size",
    y = "Computation time (minutes)",
    color = "Growth model - Dupl. rate"
  ) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Border - 0.1",
               "Border - 0.2",
               "Border - 0.4",
               "Homog. - 0.1",
               "Homog. - 0.2",
               "Homog. - 0.4")
  )

save_image(section, "size-vs-comp_time_rate", extension=img_ext,
           directory=image_dir)

cat("Alive cells-vs-time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Dup.~rate&$1/N$&$N$&$N\\log{N}$&",
           "$N^{3/2}$&$N^{2}$&$e^N$\\\\\n",
           "\\midrule\n"))
for (growth_rate_v in growth_rates) {
  for (growth_model in c(TRUE,FALSE)) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border-driven"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(total_result %>%
                                    filter(growth_rate==growth_rate_v &
                                           border_growth_model == growth_model),
                                  "num_of_alive_cells", "elapsed_time")
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name, growth_rate_v))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")
