if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

devtools::load_all()
library(peakRAM)
library(dplyr)
library(ggplot2)
library(scales)
library(latex2exp)
library(forcats)
library(cowplot)
library(fs)
library(GenomicAlignments)
library(Rsamtools)
library(parallel)

source("global.R")

get_total_genome_coverage <- function(sam_dir, n_cores = detectCores() - 1)
{
  
  sam_files <- list.files(sam_dir, pattern = "\\.sam$", full.names = TRUE)
  if (length(sam_files) == 0) stop("No SAM files found.")
  
  bam_files <- mclapply(sam_files, function(f) {
    out_bam <- sub("\\.sam$", ".sorted.bam", f)
    
    if(!file.exists(out_bam)) {
      system(paste("samtools sort -@ 1 -o", out_bam, f))
      system(paste("samtools index", out_bam))
    }
    return(out_bam)
  }, mc.cores = n_cores) |> unlist()
  
  stats_list <- lapply(bam_files, function(b) {
    aln <- readGAlignments(b)
    cov <- coverage(aln)
    
    sum_depth <- sum(as.numeric(unlist(cov)))
    total_width <- sum(as.numeric(runLength(unlist(cov))))
    
    return(c(sum_depth = sum_depth, width = total_width))
  })
  
  final_matrix <- do.call(rbind, stats_list)
  total_sum_depth <- sum(final_matrix[, "sum_depth"])
  total_genome_width <- sum(final_matrix[, "width"])
  
  final_coverage <- total_sum_depth / total_genome_width
  
  return(final_coverage)
}

output_file <- file.path(result_dir, "sequencing-time.rds")

setwd(working_dir)

sample_forest_dir <- "sample_forests"

section <- "sequencing"
image_ext <- "png"

num_of_samples <- 10
measure_every <- 1
step_size <- 400000
num_of_steps <- 5
num_of_repetitions <- 5

AICc_delta_threshold <- 1

mutation_rates <- c(5e-10, 7.5e-10, 1.0e-9, 1.25e-9, 1.5e-9)
CNAs <- seq(0, 1000, by = 200)
CNA_size <- 10000
CNA_pos <- 170*60
border_growth_list <- c(TRUE, FALSE)
sample_on_edge_list <- c(TRUE, FALSE)
insert_size_list <- c(0, 1000, 2000, 3000, 4000, 5000)
read_size_list <- c(50, 100, 150, 200, 250)
write_SAM_list <- c(TRUE, FALSE)
default_coverage <- 0.1
coverages <- seq(0.2, 1.0, by=0.2)

AICc_delta_threshold <- 1

#border_growth_list <- c(TRUE)
#sample_on_edge_list <- c(TRUE)

num_of_measures <- num_of_samples/measure_every

seq_result <- NULL
if (file.exists(output_file)) {
  seq_result <- readRDS(output_file)
}

test_df <- rbind(expand.grid(sample_on_edge_list, num_of_steps,
                       	     1:num_of_repetitions, border_growth_list,
                       	     mutation_rates, 0, num_of_measures,
                       	     0.4, 0.01, 1000, 1000, default_coverage,
                             TRUE, 0, 150),
		             expand.grid(sample_on_edge_list, num_of_steps,
		                         1:num_of_repetitions, border_growth_list,
		                         1e-9, 0, c(2,4,6,8,10),
		                         0.4, 0.01, 1000, 1000, default_coverage,
		                         TRUE, 0, 150),
		             expand.grid(sample_on_edge_list, num_of_steps,
		                         1:num_of_repetitions, border_growth_list,
		                         1e-9, 0, num_of_measures,
		                         0.4, 0.01, 1000, 1000, default_coverage,
		                         write_SAM_list, 0, 150),
		             expand.grid(sample_on_edge_list, num_of_steps,
		                         1:num_of_repetitions, border_growth_list,
		                         1e-9, 0, num_of_measures,
		                         0.4, 0.01, 1000, 1000, default_coverage,
		                         TRUE, 0, read_size_list),
		             expand.grid(sample_on_edge_list, num_of_steps,
		                         1:num_of_repetitions, border_growth_list,
		                         1e-9, 0, num_of_measures,
		                         0.4, 0.01, 1000, 1000, default_coverage,
		                         TRUE, insert_size_list, 150),
            		 expand.grid(sample_on_edge_list, num_of_steps,
                             1:num_of_repetitions, border_growth_list,
            		             1e-9, 0, num_of_measures,
            		             0.4, 0.01, 1000, 1000, coverages,
            		             TRUE, 0, 150))

m_engine <- MutationEngine(setup_code = "GRCh38",
                           context_sampling=20,
                           max_repetition_storage=2000000)

genome_info <- m_engine$get_genome_info() %>%
  mutate(allelic_size = size*num_of_alleles)

DNA_amount <- sum(genome_info$allelic_size)

for (num_of_CNAs in CNAs) {
  rate <- 1.0e-9*DNA_amount/(DNA_amount+CNA_size*as.double(num_of_CNAs))
  test_df <- rbind(test_df,
                   expand.grid(sample_on_edge_list, num_of_steps,
                               1:num_of_repetitions, border_growth_list,
                               rate, num_of_CNAs, num_of_measures,
                               0.4, 0.01, 0, 0, default_coverage,
                               TRUE, 0, 150))
}

test_df <- test_df %>% distinct()

names(test_df) <- c("sample_on_edge", "num_of_steps", "repetition_id",
                    "border_growth", "rate", "num_of_CNAs", "measure_id",
                    "duplication_rate", "death_rate", "pre_neo_SNV",
                    "pre_neo_ID", "coverage", "write_SAM", "insert_size",
                    "read_size")

set.seed(0)

test_order <- sample(1:nrow(test_df))

#test_order <- list()

saving <- TRUE

current_test <- 0
for (test_idx in test_order) {
  current_test <- current_test + 1
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
  coverage_v <- test_df[test_idx, 12]
  write_SAM_v <- test_df[test_idx, 13]
  insert_size_v <- test_df[test_idx, 14]
  read_size_v <- test_df[test_idx, 15]
  
  tumour_size <- i*step_size
  
  samples_in_forest <- k*measure_every
  
  to_measure <- is.null(seq_result)
  
  if (!to_measure) {
    in_df <- nrow(seq_result %>%
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
                             pre_neo_ID == pre_neo_ID_v & 
                             aimed_coverage == coverage_v &
                             write_SAM == write_SAM_v &
                             insert_size == insert_size_v &
                             read_size == read_size_v))
    
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
    cat(paste("  coverage:", coverage_v, "\n"))
    cat(paste("  write_SAM:", write_SAM_v, "\n"))
    cat(paste("  insert_size:", insert_size_v, "\n"))
    cat(paste("  read_size:", read_size_v, "\n"))

    set.seed(i*k)
    
    m_engine <- MutationEngine(setup_code = "GRCh38",
                               context_sampling=20,
                               max_repetition_storage=2000000)
    
    driver_list <- list()
    if (num_of_CNAs_v>0) {
      for (l in 1:num_of_CNAs_v) {
        driver_list <- append(driver_list,
                              CNA(type = "A", chr = "1", chr_pos = CNA_pos,
                                  len = CNA_size))
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

    pf <- m_engine$place_mutations(forest, pre_neo_SNV_v, pre_neo_ID_v)
    nodes <- pf$get_nodes()
    if (nrow(nodes) != nrow(forest$get_nodes())) {
      stop("Diverso numero di nodi")
    }
    
    sample_info <- forest$get_samples_info()
    sample_name <- sample_info[nrow(sample_info), "name"]
    pf <- pf$get_subforest_for(sample_name)

    nodes <- pf$get_nodes()
    mut_stats <- pf$get_mutation_statistics()
    sample_info <- pf$get_samples_info()
    
    sample_nodes <- nodes %>% filter(!is.na(sample))
    depth_stats <- sample_nodes %>%
      summarize(mean_depth=mean(depth), stdev_depth=sd(depth),
                height=max(depth), simulation_time=max(birth_time))
    
    mut_stats <- pf$get_mutation_statistics()
    
    total_SID_stats <- inner_join(sample_nodes, mut_stats,
                                  by="cell_id") %>%
      summarise(total=sum(total_SIDs),
                mean=mean(total_SIDs),
                stdev=sd(total_SIDs))

    basic_seq <- BasicIlluminaSequencer(4e-3)

    seq_dir <- "seq_test"
    
    if (dir.exists(seq_dir)) {
      unlink(seq_dir, recursive=TRUE)
    }
    
    measure <- peakRAM(seq_output <- simulate_seq(pf, sequencer = basic_seq,
                                                  coverage = coverage_v,
                                                  write_SAM = write_SAM_v,
                                                  read_size = read_size_v,
                                                  insert_size_mean = insert_size_v,
                                                  with_normal_sample = FALSE,
                                                  output_dir = seq_dir))

    if (write_SAM_v) {
      dir_data <- dir_info(seq_dir, recurse = TRUE)
      SAM_data_v <- sum(dir_data$size)
      real_coverage <- get_total_genome_coverage(seq_dir)
    } else {
      SAM_data_v <- 0
      real_coverage <- NA
    }
    
    if (dir.exists(seq_dir)) {
      unlink(seq_dir, recursive=TRUE)
    }
    
    result <- data.frame(elapsed_time=measure$Elapsed_Time_sec,
                         tumour_size = i*step_size,
                         repetition_id = j,
                         num_of_samples = samples_in_forest,
                         seq_samples = 1,
                         mutation_rate = rate,
                         num_of_CNAs = num_of_CNAs_v,
                         border_growth = border_growth_v,
                         sampled_on_edge = sample_on_edge,
                         duplication_rate = dup_rate_v,
                         death_rate = death_rate_v,
                         pre_neo_SNV = pre_neo_SNV_v,
                         pre_neo_ID = pre_neo_ID_v,
                         aimed_coverage = coverage_v,
                         write_SAM = write_SAM_v,
                         read_size = read_size_v,
                         insert_size = insert_size_v,
                         eval_coverage = real_coverage,
                         SAM_size = SAM_data_v,
                         num_of_nodes = nrow(nodes),
                         num_of_leaves = sum(sample_info$tumour_cells), 
                         forest_height = depth_stats$height,
                         total_somatic = total_SID_stats$total,
                         total_somatic_mean = total_SID_stats$mean,
                         total_somatic_stdev = total_SID_stats$stdev,
                         sample_depth_mean = depth_stats$mean,
                         sample_depth_stdev = depth_stats$stdev)
    
    if (saving) {
      if (is.null(seq_result)) {
        seq_result <- result
      } else {
        seq_result <- rbind(seq_result, result)
      }
      saveRDS(seq_result, output_file)
    }
  }
}

plot_data <- seq_result %>% filter(num_of_CNAs==0 &
                                     num_of_samples == 10 &
                                     tumour_size == 2e+06 &
                                     pre_neo_SNV == 1000 &
                                     pre_neo_ID == 1000 &
                                     aimed_coverage == 0.1 &
                                     insert_size == 0 &
                                     read_size == 150 &
                                     write_SAM == TRUE)

plot1 <- ggplot(data = plot_data,
                mapping = aes(x = mutation_rate, y = elapsed_time,
                              color = interaction(border_growth,
                                                  sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2e-9) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Mutation rate",
    y = "Sequencing time (seconds)"
  ) +
  scale_x_continuous(breaks = mutation_rates) +
  geom_smooth(method = "lm", formula = y ~ x, se=TRUE) +
  theme_minimal() +
  scale_color_discrete(
    name = NULL, 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  ) +
  theme(legend.position = "bottom")

cat("\nMutation rates-vs-Time\n")

cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/r_{m}$&$r_{m}$&$r_{m}\\log{r_{m}}$",
           "&$r_{m}^{3/2}$&$r_{m}^{2}$&$e^{r_{m}}$\\\\\n",
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

  results <- compare_complexity(plot_data %>%
                                  filter(border_growth == growth_model &
                                         sampled_on_edge == sampled_on_edge_v),
                                "mutation_rate", "elapsed_time")
  
  latex_AICc_row(results$AICc, AICc_delta_threshold,
                 first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

plot_data <- seq_result %>% filter(num_of_CNAs==0 &
                                     tumour_size == 2e+06 &
                                     num_of_samples == 10 &
                                     pre_neo_SNV == 1000 &
                                     pre_neo_ID == 1000 &
                                     aimed_coverage > 0.1 &
                                     mutation_rate == 1e-9 &
                                     insert_size == 0 &
                                     read_size == 150 &
                                     write_SAM == TRUE)

plot2 <- ggplot(data = plot_data,
                mapping = aes(x = aimed_coverage, y = elapsed_time,
                              color = interaction(border_growth,
                                                  sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 0.15) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Coverage",
    y = "Sequencing time (seconds)"
  ) +
  geom_smooth(method = "lm", formula = y ~ x, se=TRUE) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

cat("\nCoverage-vs-Time\n")
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
    
    results <- compare_complexity(plot_data %>%
                                    filter(border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "aimed_coverage", "elapsed_time")
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

plot_data <- seq_result %>% filter(tumour_size == 2e+06 &
                                     num_of_samples == 10 &
                                     pre_neo_SNV == 1000 &
                                     pre_neo_ID == 1000 &
                                     aimed_coverage == 0.1 &
                                     insert_size == 0 &
                                     num_of_CNAs == 0 &
                                     write_SAM == TRUE &
                                     mutation_rate == 1e-9)

plot3 <- ggplot(data = plot_data,
                mapping = aes(x = read_size,
                              y = elapsed_time,
                              color = interaction(border_growth,
                                                  sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 10) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Read size",
    y = "Sequencing time (seconds)"
  ) +
  geom_smooth(method = "lm", formula = y ~ I(1/x), se=TRUE) +
  #scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

cat("\nRead size-vs-Time\n")
cat(paste0("\\begin{tabular}{@{}llllllllll@{}}\n\\toprule\n",
           "Growth model&Sampling Pos.&$1/\\ell_{r}$&$\\ell_{r}$&$\\ell_{r}\\log{\\ell_{r}}$",
           "&$\\ell_{r}^{3/2}$&$\\ell_{r}^{2}$&$e^{\\ell_{r}}$\\\\\n",
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
    
    results <- compare_complexity(plot_data %>%
                                    filter(border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "read_size", "elapsed_time")
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

plot_data <- seq_result %>%
  filter(tumour_size == 2e+06 & num_of_samples == 10 &
           pre_neo_SNV == 0 & pre_neo_ID == 0 &
           aimed_coverage == 0.1 &
           insert_size == 0 & read_size == 150 &
           num_of_CNAs %in% CNAs &
           write_SAM == TRUE) %>%
  mutate(DNA_quantity = sum(genome_info$num_of_alleles) + num_of_CNAs)

plot4 <- ggplot(data = plot_data,
                mapping = aes(x = DNA_quantity, y = elapsed_time,
                              color = interaction(border_growth,
                                                  sampled_on_edge))) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 20) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Number of driver CNAs",
    y = "Sequencing time (seconds)"
  ) +
  geom_smooth(method = "lm", formula = y ~ x, se=TRUE) +
  theme_minimal() +
  scale_color_discrete(
    name = "Growth - Position", 
    labels = c("Homog. - Centre",
               "Border - Centre",
               "Homog. - Edge",
               "Border - Edge")
  )

cat("\nNum. of CNA-vs-Time\n")

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
    
    results <- compare_complexity(plot_data %>%
                                    filter(border_growth == growth_model &
                                             sampled_on_edge == sampled_on_edge_v),
                                  "num_of_CNAs", "elapsed_time")
    
    latex_AICc_row(results$AICc, AICc_delta_threshold,
                   first_columns=c(growth_model_name,sampled_on_edge_name))
  }
}
cat("\\bottomrule\n\\end{tabular}\n")

td_space = 10
right_space = 15
left_space = 15

shared_legend <- get_legend(
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot1 <- plot1 + theme(legend.position = "none",
                       plot.margin = margin(t = td_space, r = 0,
                                            b = td_space, l = left_space))
plot2 <- plot2 + theme(legend.position = "none", axis.title.y = element_blank(),
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))
plot3 <- plot3 + theme(legend.position = "none",
                       plot.margin = margin(t = td_space, r = 0,
                                            b = td_space, l = left_space))
plot4 <- plot4 + theme(legend.position = "none", axis.title.y = element_blank(),
                       plot.margin = margin(t = td_space, r = right_space,
                                            b = td_space, l = left_space))

final_plot <- plot_grid(plot1, plot2, plot3, plot4,
                       labels = c('A', 'B', 'C', 'D'), ncol = 2,
                       label_size = 12)

final_plot <- plot_grid(final_plot, shared_legend, ncol = 1,
                        rel_heights = c(1, .15))

save_image(section, "overall", width=12, height=8, extension=image_ext,
           directory=image_dir)

cat("\n")


