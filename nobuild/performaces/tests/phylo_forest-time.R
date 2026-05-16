devtools::load_all()
library(peakRAM)
library(dplyr)
library(ggplot2)
library(scales)
library(latex2exp)
library(forcats)
library(cowplot)

source("global.R")

if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

output_file <- file.path(result_dir, "phylo_forests-time.rds")

setwd(working_dir)

sample_forest_dir <- "sample_forests"

section <- "phylo"
image_ext <- "png"

num_of_samples <- 10
measure_every <- 1
step_size <- 400000
num_of_steps <- 5
num_of_repetitions <- 5

mutation_rates <- c(5e-10, 7.5e-10, 1.0e-9, 1.25e-9, 1.5e-9)
CNAs <- 0:10
border_growth_list <- c(TRUE, FALSE)
sample_on_edge_list <- c(TRUE, FALSE)

AICc_delta_threshold <- 1

num_of_measures <- num_of_samples/measure_every

forest_result <- NULL

if (file.exists(output_file)) {
  forest_result <- readRDS(output_file)
}

test_df <- rbind(expand.grid(sample_on_edge_list, num_of_steps,
                             1:num_of_repetitions, border_growth_list,
                             1.0e-9, CNAs, num_of_measures,
                             0.4, 0.01, 0, 0),
		             expand.grid(TRUE, num_of_steps,
                             1:num_of_repetitions, border_growth_list,
                             1.0e-9, 0, num_of_measures,
                             0.4, 0.01, 1000, 1000),
                 expand.grid(TRUE, 1:num_of_steps,
                             1:num_of_repetitions, TRUE,
                             1.0e-9, 0, num_of_measures,
                             0.4, 0.01, 1000, 1000),
                 expand.grid(sample_on_edge_list, num_of_steps,
                             1:num_of_repetitions, border_growth_list,
                             mutation_rates, 0, num_of_measures,
                             0.4, 0.01, 1000, 1000),
                expand.grid(sample_on_edge_list, num_of_steps,
                             1:num_of_repetitions, border_growth_list,
                             1.0e-9, 0, 1:num_of_measures,
                             0.4, 0.01, 1000, 1000)
                ) %>%
  distinct()


forest_filename <- get_sample_forest_filename(num_of_steps, 1, 1,
                                              TRUE, TRUE, 0.4, 0.01)

forest_file <- file.path(sample_forest_dir, forest_filename)

forest <- load_sample_forest(forest_file)

depth_stats <- (forest$get_nodes() %>% filter(!is.na(sample)) %>%
                  summarize(mean_depth=mean(depth), sum_height=sum(depth)))
sum_height_model <- depth_stats$sum_height

test_df <- NULL

cat("Setting up tests")

for (j in 1:num_of_repetitions) {
  for (sample_on_edge_v in sample_on_edge_list) {
    cat(".")
    for (border_growth_v in border_growth_list) {
      for (k in 1:num_of_measures) {
        forest_filename <- get_sample_forest_filename(num_of_steps, k, j,
                                                      border_growth_v,
                                                      sample_on_edge_v, 0.4,
                                                      0.01)
        
        forest_file <- file.path(sample_forest_dir, forest_filename)
        
        forest <- load_sample_forest(forest_file)

        depth_stats <- (forest$get_nodes() %>% filter(!is.na(sample)) %>%
          summarize(mean_depth=mean(depth), sum_height=sum(depth)))
        sum_height <- depth_stats$sum_height
        
        mutation_rate <- 1.0e-9*(sum_height_model/sum_height)
  
        test_df <- rbind(test_df,
                         expand.grid(sample_on_edge_v, num_of_steps,
                                     j, border_growth_v, mutation_rate,
                                     0, k, 0.4, 0.01, 0, 0))  %>%
          distinct()
      }
    }
  }
}

cat("done\n")

for (i in 1:7) {
  test_df <- rbind(test_df,
                   expand.grid(sample_on_edge_list, num_of_steps,
                               1:num_of_repetitions, border_growth_list,
                               1.0e-9, 0, num_of_measures,
                               0.4, 0.01, i*1000, i*1000))  %>%
    distinct()
}

names(test_df) <- c("sample_on_edge", "num_of_steps", "repetition_id",
                    "border_growth", "rate", "num_of_CNAs", "measure_id",
                    "duplication_rate", "death_rate", "pre_neo_SNV",
                    "pre_neo_ID")

set.seed(0)

test_order <- sample(1:nrow(test_df))

warmup <- TRUE
saving <- TRUE

if (warmup) {
  warmup_counter <- 3
} else {
  warmup_counter <- 1
}

current_test <- 1
for (test_idx in test_order) {
  sample_on_edge <- test_df[test_idx, 1]
  i <- test_df[test_idx, 2]
  j <- test_df[test_idx, 3]
  border_growth_v <- test_df[test_idx, 4]
  rate <- test_df[test_idx, 5]
  num_of_CNAs_v <- test_df[test_idx, 6]
  k <- test_df[test_idx, 7]
  dup_rate_v <- test_df[test_idx, 8]
  death_rate_v <- test_df[test_idx, 9]
  pre_neo_SNV_v <- test_df[test_idx, 10]
  pre_neo_ID_v <- test_df[test_idx, 11]
        
  tumour_size <- i*step_size

  samples_in_forest <- k*measure_every
  
  to_measure <- is.null(forest_result)
  
  if (!to_measure) {
    in_df <- nrow(forest_result %>%
                        filter(tumour_size == i*step_size &
                                 repetition_id == j &
                                 num_of_samples == samples_in_forest &
                                 mutation_rate == rate &
                                 num_of_CNAs == num_of_CNAs_v &
                                 border_growth == border_growth_v &
                                 sampled_on_edge == sample_on_edge &
                                 duplication_rate == dup_rate_v &
                                 death_rate == death_rate_v &
                                 pre_neo_SNV == pre_neo_SNV_v &
                                 pre_neo_ID == pre_neo_ID_v))
    
    to_measure <- (in_df == 0)
  }
  
  if (to_measure) {
    cat(paste0("Setup #", current_test," of ", nrow(test_df), ":\n"))
    cat(paste("  samples_in_forest:", samples_in_forest, "\n"))
    cat(paste("  tumour_size:", tumour_size, "\n"))
    cat(paste("  repetition_id:", j, "\n"))
    cat(paste("  rate:", rate, "\n"))
    cat(paste("  num. of CNAs:", num_of_CNAs_v, "\n"))
    cat(paste("  border_growth:", border_growth_v, "\n"))
    cat(paste("  sample_on_edge:", sample_on_edge, "\n"))
    cat(paste("  duplication_rate:", dup_rate_v, "\n"))
    cat(paste("  death_rate:", death_rate_v, "\n"))
    cat(paste("  pre_neo_SNV:", pre_neo_SNV_v, "\n"))
    cat(paste("  pre_neo_ID:", pre_neo_ID_v, "\n"))

    while (warmup_counter>0) {
      set.seed(i*k)

      m_engine <- MutationEngine(setup_code = "GRCh38",
                                 context_sampling=20,
                                 max_repetition_storage=2000000)

      chr_1_info <- m_engine$get_genome_info() %>% filter(chr == "1")
      driver_list <- list()
      if (num_of_CNAs_v>0) {
        for (l in 1:num_of_CNAs_v) {
          driver_list <- append(driver_list,
                                CNA(type = "A",chr = chr_1_info$chr,
                                    chr_pos = 1, len = chr_1_info$size))
        }
      }
      
      m_engine$add_mutant(mutant_name = "Clone 1",
                          passenger_rates = list(SNV = rate,
                                                 indel = rate),
                          drivers = driver_list)
      
      m_engine$add_exposure(c(ID2 = 0.6, ID13 = 0.2, ID21 = 0.2,
                              SBS13 = 0.2, SBS1 = 0.8))
  
      forest_filename <- get_sample_forest_filename(i, k, j, border_growth_v,
                                                    sample_on_edge, dup_rate_v,
                                                    death_rate_v)
  
      forest_file <- file.path(sample_forest_dir, forest_filename)
      
      forest <- load_sample_forest(forest_file)
      
      measure <- peakRAM(pf <- m_engine$place_mutations(forest, pre_neo_SNV_v,
                                                        pre_neo_ID_v))
      
      nodes <- pf$get_nodes()
      if (nrow(nodes) != nrow(forest$get_nodes())) {
        stop("Diverso numero di nodi")
      }
      mut_stats <- pf$get_mutation_statistics()
  
      sample_nodes <- nodes %>% filter(!is.na(sample))
      depth_stats <- sample_nodes %>%
        summarize(mean_depth=mean(depth), stdev_depth=sd(depth),
                  height=max(depth), simulation_time=max(birth_time))
  
      new_SID_stats <- mut_stats %>%
        summarise(mean=mean(new_SIDs),
                  stdev=sd(new_SIDs),
                  total=sum(new_SIDs))
  
      total_SID_stats <- inner_join(sample_nodes, mut_stats,
                                   by="cell_id") %>%
        summarise(total=sum(total_SIDs),
                  mean=mean(total_SIDs),
                  stdev=sd(total_SIDs))
  
      sample_info <- forest$get_samples_info()
  
      result <- data.frame(elapsed_time=measure$Elapsed_Time_sec,
                           repetition_id=j,
                           tumour_size = i*step_size,
                           num_of_samples = nrow(sample_info),
                           num_of_nodes = nrow(nodes),
                           num_of_leaves = sum(sample_info$tumour_cells),
                           total_somatic = total_SID_stats$total,
                           total_somatic_mean = total_SID_stats$mean,
                           total_somatic_stdev = total_SID_stats$stdev,
                           distinct_somatic = new_SID_stats$total,
                           distinct_somatic_mean = new_SID_stats$mean,
                           distinct_somatic_stdev = new_SID_stats$stdev, 
                           forest_height = depth_stats$height,
                           sample_depth_mean = depth_stats$mean,
                           sample_depth_stddev = depth_stats$stdev,
                           simulation_time = depth_stats$simulation_time,
                           mutation_rate = rate,
                           num_of_CNAs = num_of_CNAs_v,
                           border_growth = border_growth_v,
                           sampled_on_edge = sample_on_edge,
                           duplication_rate = dup_rate_v,
                           death_rate = death_rate_v,
                           pre_neo_SNV = pre_neo_SNV_v,
                           pre_neo_ID = pre_neo_ID_v)
      
      rm(pf)
      rm(forest)
      rm(nodes)
      rm(sample_nodes)
      
      warmup_counter <- warmup_counter-1
    }

    warmup_counter <- 1

    if (saving) {
      if (is.null(forest_result)) {
        forest_result <- result
      } else {
        forest_result <- rbind(forest_result, result)
      }
      saveRDS(forest_result, output_file)
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

plot1 <- ggplot(data = forest_result %>% filter(num_of_CNAs==0 &
                                                   tumour_size == 2e+06 &
                                                    pre_neo_SNV == 0 &
                                                    pre_neo_ID == 0 &
                                                  (mutation_rate != 1e-9 |
                                                     (border_growth == TRUE &
                                                        sampled_on_edge == TRUE &
                                                        num_of_samples == 1))),
       mapping = aes(x = num_of_leaves*mutation_rate, y = elapsed_time,
                     color = interaction(border_growth, sampled_on_edge))) + 
  geom_point(size = 3) + 
  labs(
    x = "Sampled cells",
    y = "Building time (seconds)"
  ) + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  #geom_smooth(method = "lm", formula = y ~ x:log(x)+x) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )


#save_image(section, "sampled_cells-vs-time", extension=image_ext,
#           directory=image_dir)

cat("\nSampled cells-vs-Time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/S$&$S$&$S\\log{S}$",
           "&$S^{3/2}$&$S^{2}$&$e^{S}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "On Edge"
  } else {
    sampled_on_edge_name <- "In Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(forest_result %>% 
                                    filter(num_of_CNAs==0 &
                                             tumour_size == 2e+06 &
                                             pre_neo_SNV == 0 &
                                             pre_neo_ID == 0  &
                                             border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "num_of_leaves", "elapsed_time",
                                  #paste0("Nodes in forest-vs-Time (",
                                  #       growth_model_name, " - ",
                                  #       sampled_on_edge_name, ")")
                                  )
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

plot1 <- ggplot(data = forest_result %>% filter(num_of_CNAs==0 &
                                                  mutation_rate == 1e-9 &
                                                  tumour_size == 2e+06 &
                                                  pre_neo_SNV == 1000 &
                                                  pre_neo_ID == 1000),
                mapping = aes(x = num_of_leaves, y = elapsed_time,
                              color = interaction(border_growth, sampled_on_edge))) + 
  geom_point(size = 3) + 
  labs(
    x = "Sampled cells",
    y = "Building time (seconds)"
  ) + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  geom_smooth(method = "lm", formula = y ~ x:log(x)+x) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )


#save_image(section, "sampled_cells-vs-time", extension=image_ext,
#           directory=image_dir)

cat("\nSampled cells-vs-Time 2\n")
cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/S$&$S$&$S\\log{S}$",
           "&$S^{3/2}$&$S^{2}$&$e^{S}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "On Edge"
  } else {
    sampled_on_edge_name <- "In Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(forest_result %>% 
                                    filter(num_of_CNAs==0 &
                                             tumour_size == 2e+06 &
                                             pre_neo_SNV == 1000 &
                                             pre_neo_ID == 1000  &
                                             border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "num_of_leaves", "elapsed_time",
                                  #paste0("Nodes in forest-vs-Time (",
                                  #       growth_model_name, " - ",
                                  #       sampled_on_edge_name, ")")
    )
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

sample_categories <- as.character(rev(unique(forest_result$num_of_samples)))

rate_categories <- as.character(sort(unique(forest_result$mutation_rate)))

forest_result$rate_categories <- as.character(forest_result$mutation_rate)

forest_result$rate_categories <- factor(forest_result$rate_categories,
                                       levels = rate_categories)

plot2 <- ggplot(forest_result %>% filter(num_of_CNAs==0 &
                                             tumour_size == 2e+6 & 
                                             num_of_samples == 10 &
                                           pre_neo_SNV == 1000 &
                                           pre_neo_ID == 1000),
                  aes(x = rate_categories, y = elapsed_time,
                      color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
               geom = "errorbar", width = 1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Mutation Rates",
    y = "Building time (seconds)"
  ) +
  #scale_x_continuous(labels = label_number(scale = 1e-9, suffix = "B")) +
  #geom_smooth(method = "lm", formula = y ~ I(x * log(x)), color = "blue") +
  geom_smooth(method = "lm", formula = y ~ x:log(x)+x, se=TRUE) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

plot2 <- ggplot(forest_result %>% filter(num_of_CNAs==0 &
                                           tumour_size == 2e+6 & 
                                           num_of_samples == 10 &
                                           pre_neo_SNV == 1000 &
                                           pre_neo_ID == 1000),
                aes(x = mutation_rate, y = elapsed_time,
                    color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 1.5e-10) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Mutation rates",
    y = "Building time (seconds)"
  ) +
  #scale_x_continuous(labels = label_number(scale = 1e-9, suffix = "B")) +
  #geom_smooth(method = "lm", formula = y ~ I(x * log(x)), color = "blue") +
  geom_smooth(method = "lm", formula = y ~ x:log(x)+x, se=TRUE) +
  theme_minimal() +
  scale_x_continuous(labels = label_scientific(digits=2)) +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

#save_image(section, "mutation_rate-vs-time", extension=image_ext,
#           directory=image_dir)

cat("\nMutation rates-vs-Time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/S$&$S$&$S\\log{S}$",
           "&$S^{3/2}$&$S^{2}$&$e^{S}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "On Edge"
  } else {
    sampled_on_edge_name <- "In Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(forest_result %>% 
                                    filter(num_of_CNAs==0 &
                                             tumour_size == 2e+6 & 
                                             num_of_samples == 10 &
                                             pre_neo_SNV == 1000 &
                                             pre_neo_ID == 1000 & 
                                             border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "mutation_rate", "elapsed_time",
                                  #paste0("Mutation rates-vs-Time (",
                                  #       growth_model_name, " - ",
                                  #       sampled_on_edge_name, ")")
                                  )

    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

preneo_data <- forest_result %>% filter(num_of_CNAs==0 &
                                          tumour_size == 2e+6 & 
                                          num_of_samples == 10 &
                                          mutation_rate == 1e-9 &
                                          (pre_neo_SNV+pre_neo_ID) %% 2000 == 0 &
                                          (pre_neo_SNV+pre_neo_ID) > 5000)

plot3 <- ggplot(preneo_data,
                aes(x = pre_neo_SNV+pre_neo_ID, y = elapsed_time,
                    color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 750) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Pre-neoplastic mutations",
    y = "Building time (seconds)"
  ) +
  geom_smooth(data = subset(preneo_data, border_growth == TRUE & 
                              sampled_on_edge == TRUE), 
              method = "lm", formula = y ~ log(x)+x) +
  geom_smooth(data = subset(preneo_data, border_growth == TRUE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ log(x)+x) +
  geom_smooth(data = subset(preneo_data, border_growth == FALSE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ log(x)+x) +
  geom_smooth(data = subset(preneo_data, border_growth == FALSE & 
                              sampled_on_edge == TRUE), 
              method = "glm", formula = y ~ log(x)+x) +
  #scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

#save_image(section, "pre-neoplastic-vs-time", extension=image_ext,
#           directory=image_dir)

cat("\nPre-neoplastic-vs-Time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/p$&$p$&$p\\log{p}$",
           "&$p^{3/2}$&$p^{2}$&$e^{p}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "On Edge"
  } else {
    sampled_on_edge_name <- "In Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(preneo_data %>% 
                                    filter(border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "pre_neo_ID", "elapsed_time",
                                  #paste0("Pre-neoplastic-vs-Time (",
                                  #       growth_model_name, " - ",
                                  #       sampled_on_edge_name, ")")
                                  )
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

CNA_data <- forest_result %>% filter(tumour_size == 2e+6 &
                                       mutation_rate == 1e-9 &
                                       pre_neo_SNV == 0 &
                                       pre_neo_ID == 0 &
                                       num_of_samples == 10 &
                                       num_of_CNAs %% 2 == 0)

label_breaks <- seq(0, 10, by=2)

plot4 <- ggplot(CNA_data,
                aes(x = num_of_CNAs, y = elapsed_time,
                    color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Number of driver CNAs",
    y = "Building time (seconds)"
  ) +
  #scale_x_continuous(labels = label_number(scale = 1e-9, suffix = "B")) +
  #geom_smooth(method = "lm", formula = y ~ I(x * log(x)), color = "blue") +
  geom_smooth(data = subset(CNA_data, border_growth == TRUE & 
                              sampled_on_edge == TRUE), 
              method = "lm", formula = y ~ x) +
  geom_smooth(data = subset(CNA_data, border_growth == TRUE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ x:log(x+1)+x) +
  geom_smooth(data = subset(CNA_data, border_growth == FALSE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ x:log(x+1)+x) +
  geom_smooth(data = subset(CNA_data, border_growth == FALSE & 
                              sampled_on_edge == TRUE), 
              method = "glm", formula = y ~ x:log(x+1)+x) +
  scale_x_continuous(breaks = label_breaks,
                     labels = as.character(label_breaks)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

#save_image(section, "CNA-vs-time", extension=image_ext,
#           directory=image_dir)

cat("\nCNA-vs-Time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/C$&$C$&$C\\log{C}$",
           "&$C^{3/2}$&$C^{2}$&$e^{C}$\\\\\n",
           "\\midrule\n"))
for (sampled_on_edge_v in sample_on_edge_list) {
  if (sampled_on_edge_v == TRUE) {
    sampled_on_edge_name <- "On Edge"
  } else {
    sampled_on_edge_name <- "In Center"
  }
  for (growth_model in border_growth_list) {
    if (growth_model == TRUE) {
      growth_model_name <- "Border"
    } else {
      
      growth_model_name <- "Homogeneous"
    }
    results <- compare_complexity(CNA_data %>% 
                                    filter(border_growth == growth_model &
                                           sampled_on_edge == sampled_on_edge_v),
                                  "num_of_leaves", "elapsed_time")
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

mut_df <- forest_result  %>% filter(num_of_CNAs==0 &
                                      tumour_size == 2e+6 &
                                      mutation_rate == 1e-9 &
                                      pre_neo_SNV == 1000 &
                                      pre_neo_ID == 1000)

plot5 <- ggplot(mut_df,
                aes(x = num_of_leaves, y = total_somatic,
                    color = interaction(border_growth, sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "errorbar", width = 750) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Sampled cells",
    y = "Somatic mutations"
  ) +
  geom_smooth(data = subset(mut_df, border_growth == TRUE & 
                              sampled_on_edge == TRUE), 
              method = "lm", formula = y ~ x) +
  geom_smooth(data = subset(mut_df, border_growth == TRUE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ x) +
  geom_smooth(data = subset(mut_df, border_growth == FALSE & 
                              sampled_on_edge == FALSE), 
              method = "glm", formula = y ~ x) +
  geom_smooth(data = subset(mut_df, border_growth == FALSE & 
                              sampled_on_edge == TRUE), 
              method = "glm", formula = y ~ x) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "K")) +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

td_space = 10
right_space = 15
left_space = 15

shared_legend <- get_legend(
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot1 <- plot1 + theme(legend.position = "none", #axis.title.y = element_blank(),
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot2 <- plot2 + theme(legend.position = "none",
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot3 <- plot3 + theme(legend.position = "none", #axis.title.y = element_blank(),
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot4 <- plot4 + theme(legend.position = "none",
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot5 <- plot5 + theme(legend.position = "none",
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))

plots_row1 <- plot_grid(plot4, plot1,
                       labels = c('A', 'B'), ncol = 2,
                       label_size = 12)
plots_row2 <- plot_grid(plot2, plot3,
                        labels = c('C', 'D'), ncol = 2,
                        label_size = 12)
plots_row3 <- plot_grid(plot5, shared_legend, ncol = 2,
                       labels = c('E',''), label_size = 12)

final_plot <- plot_grid(plot2, plot3, plot4, shared_legend, plot5, plot1,
                        labels = c('A', 'B', 'C', '', 'D', 'E'), ncol = 2,
                        label_size = 12)

save_image(section, "overall", width=12, height=8, extension=image_ext,
           directory=image_dir)

