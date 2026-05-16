devtools::load_all()
library(peakRAM)
library(dplyr)
library(ggplot2)
library(scales)
library(forcats)
library(latex2exp)
library(cowplot)

source("global.R")

get_center_sample_squares <- function(tissue_sizes, sample_size, num_of_samples)
{
  base_x <- (tissue_sizes[1]-sample_size)/2
  base_y <- (tissue_sizes[2]-sample_size)/2
  
  radius <- floor(ceiling(sqrt(num_of_samples))/2)
  
  sample_squares <- c()
  
  for (delta_x in -radius:radius) {
    for (delta_y in -radius:radius) {
      x = base_x + delta_x * (sample_size+1)
      y = base_y + delta_y * (sample_size+1)

      sample_squares <- c(sample_squares,
                          new(TissueRectangle, c(x, y),
                              c(x+sample_size, y+sample_size)))
    }
  }

  dist <- function(sample_square) {
    lc <- sample_square$lower_corner

    return(max(abs(lc[1]-base_x), abs(lc[2]-base_y)))
  }

  return(sample_squares[order(sapply(sample_squares, dist))])
}

if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

setwd(working_dir)

output_dir <- "sample_forests"

output_forest_file <- file.path(result_dir, "sample_forests-time.rds")
output_sampling_file <- file.path(result_dir, "sampling-time.rds")

num_of_samples <- 10
measure_every <- 1
step_size <- 400000
num_of_steps <- 5
sample_size <- 70
num_of_repetitions <- 5
num_of_creations <- 40
sample_on_edge_list <- c(TRUE, FALSE)
border_growth_list <- c(TRUE, FALSE)
dup_rate_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)

AICc_delta_threshold <- 1

section <- "sampling"
image_ext <- "png"

num_of_measures <- num_of_samples/measure_every

forest_result <- NULL
if (file.exists(output_forest_file)) {
  forest_result <- readRDS(output_forest_file)
}

sampling_result <- NULL
if (file.exists(output_sampling_file)) {
  sampling_result <- readRDS(output_sampling_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

sample_squares <- get_center_sample_squares(c(7000, 7000), sample_size,
                                            num_of_samples)

test_df <- rbind(expand.grid(sample_on_edge_list, c(num_of_steps),
                             1:num_of_repetitions, border_growth_list, c(0.4),
                             c(0.01)),
                 expand.grid(sample_on_edge_list, 1:num_of_steps,
                             1:num_of_repetitions, c(TRUE), c(0.4),
                             c(0.01)),
                 expand.grid(c(TRUE), c(num_of_steps),
                        1:num_of_repetitions, c(TRUE), dup_rate_list,
                        c(0.01)),
                 expand.grid(sample_on_edge_list, c(num_of_steps),
                             1:num_of_repetitions, border_growth_list,
                             dup_rate_list, c(0.01)))

for (dup_rate_v in dup_rate_list) {
  test_df <- rbind(test_df, expand.grid(c(TRUE), c(num_of_steps),
                                        1:num_of_repetitions, border_growth_list,
                                        c(dup_rate_v),
                                        c(dup_rate_v/10)))
}

test_df <- test_df %>% distinct()
names(test_df) <- c("sample_on_edge", "num_of_steps", "repetition_id",
                    "border_growth", "duplication_rate", "death_rate")

set.seed(0)

test_order <- sample(1:nrow(test_df))

test_order <- 1:nrow(test_df)

current_test <- 1
for (test_idx in test_order) {
  sample_on_edge <- test_df[test_idx, 1]
  i <- test_df[test_idx, 2]
  j <- test_df[test_idx, 3]
  border_growth_v <- test_df[test_idx, 4]
  dup_rate_v <- test_df[test_idx, 5]
  death_rate_v <- test_df[test_idx, 6]

  set.seed(j)

  tumour_size <- i*step_size

  to_measure <- is.null(forest_result)
  
  if (!to_measure) {
    recorded_measures <- nrow(forest_result %>%
                              filter(tumour_size == (i*step_size) &
                                       repetition_id == j &
                                       border_growth == border_growth_v &
                                       sampled_on_edge == sample_on_edge &
                                       duplication_rate == dup_rate_v &
                                       death_rate == death_rate_v))
    
    to_measure <- (recorded_measures < num_of_measures)
    
    if (to_measure && recorded_measures>0) {
      forest_result <- forest_result %>%
        filter(tumour_size != (i*step_size) |
                 repetition_id != j |
                 border_growth != border_growth_v |
                 sampled_on_edge != sample_on_edge |
                 duplication_rate != dup_rate_v |
                 death_rate != death_rate_v)
      
      saveRDS(forest_result, output_forest_file)
      
      sampling_result <- sampling_result %>%
        filter(tumour_size != (i*step_size) |
                 repetition_id != j |
                 border_growth != border_growth_v |
                 sampled_on_edge != sample_on_edge |
                 duplication_rate != dup_rate_v |
                 death_rate != death_rate_v)
      
      saveRDS(sampling_result, output_sampling_file)
    }
  }

  if (to_measure) {
    cat(paste0("Setup #", current_test," of ", nrow(test_df), ":\n"))
    cat(paste("  tumour_size:", tumour_size, "\n"))
    cat(paste("  repetition_id:", j, "\n"))
    cat(paste("  border_growth:", border_growth_v, "\n"))
    cat(paste("  sample_on_edge:", sample_on_edge, "\n"))
    cat(paste("  duplication_rate:", dup_rate_v, "\n"))
    cat(paste("  death_rate:", death_rate_v, "\n"))

    sim <- setup_tumour(border_growth_v, dup_rate_v, death_rate_v)
    sim <- grow_tumour_up_to(sim, tumour_size)
    
    counts <- sim$get_counts()
    
    for (k in 1:num_of_measures) {
      for (m in 1:measure_every) {
        sample_idx <- ((k-1)*measure_every+m)
        
        if (sample_on_edge) {
          bbox <- sim$search_sample(c("Clone 1" = sample_size*sample_size*0.8),
                                    sample_size, sample_size)
        } else {
          bbox <- sample_squares[[sample_idx]]
        }

        sample_name <- paste0("S_", i, "_", sample_idx, "_", j)
        measure <- peakRAM(sim$sample_cells(sample_name,
                                            bbox$lower_corner,
                                            bbox$upper_corner))

        result <- data.frame(elapsed_time=measure$Elapsed_Time_sec,
                             tumour_size=counts$counts,
                             num_of_overall_cells=counts$overall,
                             repetition_id=j,
                             border_growth = border_growth_v,
                             simulation_time=sim$get_clock(),
                             sampled_on_edge = sample_on_edge,
                             sample_index = sample_idx,
                             duplication_rate = dup_rate_v,
                             death_rate = death_rate_v
        )
        
        if (is.null(sampling_result)) {
          sampling_result <- result
        } else {
          sampling_result <- rbind(sampling_result, result)
        }
        
        saveRDS(sampling_result, output_sampling_file)
      }
      
      repeat_forest_creation <- function(sim, num_of_creations) {
        for (repeat_idx in 1:num_of_creations) {
          forest <- sim$get_sample_forest()
        }
        return(sim$get_sample_forest())
      }
      
      measure <- peakRAM(forest <- repeat_forest_creation(sim,
                                                          num_of_creations))
      
      forest_filename <- get_sample_forest_filename(i, k, j, border_growth_v,
                                                    sample_on_edge, dup_rate_v,
                                                    death_rate_v)
      forest_file <- file.path(output_dir, forest_filename)

      forest$save(forest_file)
      
      nodes <- forest$get_nodes()

      sample_info <- forest$get_samples_info()
      
      sample_nodes <- nodes %>% filter(!is.na(sample))
      depth_stats <- sample_nodes %>%
        summarize(mean_depth=mean(depth), stdev_depth=sd(depth),
                  height=max(depth), simulation_time=max(birth_time))

      result <- data.frame(elapsed_time=measure$Elapsed_Time_sec/num_of_creations,
                           tumour_size=counts$counts,
                           num_of_overall_cells=counts$overall,
                           repetition_id=j,
                           num_of_samples = nrow(sample_info),
                           num_of_nodes = nrow(nodes),
                           num_of_leaves = sum(sample_info$tumour_cells),
                           border_growth = border_growth_v,
                           forest_height = depth_stats$height,
                           sample_depth_mean = depth_stats$mean,
                           sample_depth_stddev = depth_stats$stdev,
                           simulation_time=sim$get_clock(),
                           sampled_on_edge = sample_on_edge,
                           duplication_rate = dup_rate_v,
                           death_rate = death_rate_v
                           )
  
      if (is.null(forest_result)) {
        forest_result <- result
      } else {
        forest_result <- rbind(forest_result, result)
      }
      
      saveRDS(forest_result, output_forest_file)
    }
  }
  
  current_test <- current_test + 1
}

forest_result <- forest_result %>%
  mutate(tumour_size_category = sprintf("%.2e", tumour_size))
forest_result$tumour_size_category <- fct_relevel(forest_result$tumour_size_category,
                                                  "2.00e+06", "1.60e+06",
                                                  "1.20e+06", "8.00e+05",
                                                  "4.00e+05")

forest_result$dup_rate_category <- fct_relevel(as.character(forest_result$duplication_rate),
                                           "0.1", "0.2", "0.3", "0.4", "0.5")

forest_num_of_samples <- as.character(forest_result$num_of_samples)

forest_result$num_of_samples_category <- fct_relevel(forest_num_of_samples,
                                                    "1", "2", "3", "4", "5",
                                                    "6", "7", "8", "9", "10")

my_palette <- palette_six()

plot1 <- ggplot(forest_result %>% filter(tumour_size == 2e+6 &
                                             duplication_rate == 0.4 & 
                                             death_rate == 0.01),
                  aes(x = num_of_samples_category, y = forest_height,
                      color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.75) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Number of samples",
    y = "Forest height"
  ) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

save_image(section, "samples-vs-height", extension=image_ext,
           directory=image_dir)

plot2 <- ggplot(forest_result %>% filter(tumour_size == 2e+6 &
                                             duplication_rate == 0.4 &
                                             death_rate == 0.01),
                  aes(x = num_of_nodes, y = elapsed_time,
                      color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Nodes in forest",
    y = "Building time (seconds)"
  ) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  geom_smooth(method = "lm", formula = y ~ x:log(x)+x, color = "black") +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

save_image(section, "forest_nodes-vs-time", extension=image_ext,
           directory=image_dir)

cat("Nodes-vs-time\n")

cat(paste0("\\begin{tabular}{@{}lllllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&Dup.~rate&$1/q$&$q$&$q\\log{q}$",
           "&$q^{3/2}$&$q^{2}$&$e^{q}$\\\\\n",
           "\\midrule\n"))
for (dup_rate in dup_rate_list) {
  for (sampled_on_edge_v in sample_on_edge_list) {
    if (sampled_on_edge_v == TRUE) {
      sampled_on_edge_name <- "Edge"
    } else {
      sampled_on_edge_name <- "Center"
    }
    for (growth_model in border_growth_list) {
      if (growth_model == TRUE) {
        growth_model_name <- "Border-driven"
      } else {
        
        growth_model_name <- "Homogeneous"
      }
  
      results <- compare_complexity(forest_result %>% 
                                      filter(tumour_size == 2e+6 &
                                               duplication_rate == dup_rate &
                                               death_rate == 0.01 &
                                               border_growth == growth_model & 
                                               sampled_on_edge == sampled_on_edge_v),
                                    "num_of_nodes", "elapsed_time")
      latex_AICc_row(results$AICc, AICc_delta_threshold,
                     first_columns=c(growth_model_name,sampled_on_edge_name,
                                     dup_rate))
    }
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

cumulative_sampling_result <- sampling_result %>%
  group_by(tumour_size, repetition_id, border_growth, sampled_on_edge,
           duplication_rate) %>%
  mutate(cumulative_time = sapply(sample_index, function(x){
    sum(elapsed_time[sample_index <= x])
  })
  ) %>%
  ungroup()

forest_result$num_of_samples <- as.integer(forest_result$num_of_samples)

sampling_time_df <- left_join(cumulative_sampling_result, forest_result,
                              by=c(sample_index="num_of_samples",
                                   border_growth="border_growth",
                                   sampled_on_edge="sampled_on_edge",
                                   repetition_id="repetition_id",
                                   tumour_size="tumour_size",
                                   duplication_rate="duplication_rate",
                                   death_rate="death_rate")) %>%
  select(cumulative_time, sample_index, border_growth, sampled_on_edge,
         repetition_id, tumour_size, duplication_rate, death_rate,
         num_of_leaves)

plot3 <- ggplot(sampling_time_df %>% filter(tumour_size == 2e+6 &
                                                duplication_rate == 0.4 &
                                                death_rate == 0.01),
                  aes(x = as.integer(num_of_leaves), y = cumulative_time,
                      color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Sampled cells",
    y = "Sampling time (seconds)"
  ) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1)) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

save_image(section, "sampled_cells-vs-sampling_time", extension=image_ext,
           directory=image_dir)

cat("\nSamples-vs-time\n")
cat(paste0("\\begin{tabular}{@{}lllllllllll@{}}\n",
           "\\toprule\n",
           "Growth model&Sampling Pos.&$1/n_s$&$n_s$&$n_s\\log{n_s}$",
           "&$n_{s}^{3/2}$&$n_{s}^{2}$&$e^{n_s}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "Edge"
  } else {
    sampled_on_edge_name <- "Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border-driven"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(sampling_time_df %>% 
                                    filter(tumour_size == 2e+6 &
                                             duplication_rate == 0.4 &
                                             death_rate == 0.01 &
                                             border_growth == growth_model & 
                                             sampled_on_edge==sampled_on_edge_v),
                                  "num_of_leaves", "cumulative_time")
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n\n")

new_rows <- forest_result %>% filter(tumour_size == 2e+6 &
                         border_growth == TRUE &
                         sampled_on_edge == TRUE &
                         num_of_samples == 10 &
                         death_rate == 0.01 &
                         death_rate*10 == duplication_rate)

new_rows$rate_class <- 2

rate_result <- forest_result %>% mutate(rate_class = case_when(
  death_rate*10 == duplication_rate ~ as.integer(1),
  death_rate == 0.01 ~ as.integer(2),
  TRUE ~ as.integer(0)
))

rate_result <-rbind(rate_result, new_rows)

rate_result$rate_class <- as.character(rate_result$rate_class)

rate_result$rate_class_category <- fct_relevel(rate_result$rate_class,
                                               "2", "1")

plot4 <- ggplot(rate_result %>% filter(tumour_size == 2e+6 &
                                            border_growth == TRUE &
                                            sampled_on_edge == TRUE &
                                            num_of_samples == 10 &
                                          rate_class != 0),
                 aes(x = dup_rate_category, y = num_of_nodes,
                     color = rate_class_category, group=rate_class_category)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.75) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Duplication rate",
    y = "Nodes in forest"
  ) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  theme_minimal() + 
  scale_color_discrete(
      name = "Death rate",
      labels = c(TeX(r"($0.01$)"),
                 TeX(r"($r_{\lambda}/10$)"))
  )

save_image(section, "dup_rate-vs-nodes", extension=image_ext,
           directory=image_dir)

plot5 <- ggplot(forest_result %>% filter(tumour_size == 2e+6 &
                                             duplication_rate == 0.4 &
                                             death_rate == 0.01),
                  aes(x = num_of_leaves, y = num_of_nodes,
                      color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Sampled cells",
    y = "Nodes in forest"
  ) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1)) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

save_image(section, "sampled_cells-vs-nodes", extension=image_ext,
           directory=image_dir)


td_space = 10
right_space = 15
left_space = 15

shared_legend <- get_legend(
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot1 <- plot1 + theme(plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot5 <- plot5 + theme(legend.position = "none")

plot5 <- plot_grid(plot5, shared_legend,  ncol = 2,
                   rel_widths = c(0.72, .25))
plot5 <- plot5 + theme(plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot4 <- plot4 + theme(plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot2 <- plot2 + theme(plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))

final_plot <- plot_grid(plot1, plot5, plot4, plot2,
                       labels = c('A', 'B', 'C', 'D'), ncol = 2,
                       label_size = 12)

save_image(section, "overall", width=12, height=8, extension=image_ext,
           directory=image_dir)
