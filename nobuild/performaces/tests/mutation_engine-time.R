devtools::load_all()
library(peakRAM)
library(dplyr)
library(ggplot2)
library(scales)
library(latex2exp)
library(forcats)
library(cowplot)
library(getPass)

source("global.R")

if (!(exists("result_dir") & exists("working_dir") & exists("image_dir"))) {
  stop(paste("This script is not meant to be invocated directly.",
             "Execute \"../performace_tests.\"."))
}

output_file <- file.path(result_dir, "mutation_engine-time.rds")

section <- "mutation_engine"
image_ext <- "png"

print_table <- FALSE
setup_ProCESS_code <- "GRCh38" 
num_of_repetitions <- 5
context_samplings <- c(20, 40, 60, 80, 100)
max_repetition_storages <- c(4e+5, 8e+5, 1.2e+6, 1.6e+6, 2e+6)

AICc_delta_threshold <- 1

data_result <- NULL
if (file.exists(output_file)) {
  data_result <- readRDS(output_file)
}

test_df <- rbind(expand.grid(c(FALSE, TRUE),1:num_of_repetitions,
                             context_samplings, min(max_repetition_storages)),
                 expand.grid(c(FALSE, TRUE),1:num_of_repetitions,
                             min(context_samplings),
                             max_repetition_storages)) %>%
  distinct()

names(test_df) <- c("load_saved", "repetition_id",
                    "context_sampling", "max_repetition_storage")

warmup <- FALSE
saving <- TRUE

if (warmup) {
  warmup_counter <- 3
} else {
  warmup_counter <- 1
}

setwd(working_dir)

COSMIC_account <- NULL
if (!dir.exists(setup_ProCESS_code))
{
  if (Sys.getenv("RSTUDIO") == "1") {
    COSMIC_email <- readline("Enter your COSMIC username: ")
  } else {
    cat("Enter your COSMIC username: ")
    COSMIC_email <- readLines(con = "stdin", n = 1)
  }
  COSMIC_pw <- getPass(msg = "Enter your COSMIC password: ")
  
  MutationEngine(setup_code = setup_ProCESS_code,
                COSMIC_account = list(email=COSMIC_email, password=COSMIC_pw))
}

current_test <- 1
for (test_idx in 1:nrow(test_df)) {
  l_saved <- test_df[test_idx, 1]
  r_id <- test_df[test_idx, 2]
  c_sampling <- test_df[test_idx, 3]
  max_r_storage <- test_df[test_idx, 4]

  to_measure <- is.null(data_result)
  if (!to_measure) {
    in_df <- nrow(data_result %>%
                        filter(load_saved == l_saved &
                                 repetition_id == r_id &
                                 context_sampling == c_sampling &
                                 max_repetition_storage == max_r_storage))
    
    to_measure <- (in_df == 0)
  }
  
  if (to_measure) {
    cat(paste0("Setup #", current_test," of ", nrow(test_df), ":\n"))
    cat(paste("  load_saved:", l_saved, "\n"))
    cat(paste("  repetition_id:", r_id, "\n"))
    cat(paste("  context_sampling:", c_sampling, "\n"))
    cat(paste("  max_repetition_storage:", max_r_storage, "\n"))
    
    while (warmup_counter>0) {
      
      set.seed(c_sampling*max_r_storage*r_id)
      
      cif <- file.path(setup_ProCESS_code,
                       paste0("context_index_", c_sampling, ".cif"))
      rsif <- file.path(setup_ProCESS_code,
                        paste0("rs_index_50_", as.integer(max_r_storage),
                               ".rsif"))
      dat <- file.path(setup_ProCESS_code, "germline_data",
                       "germline_NA18486.dat")
      if (l_saved) {
        if (!(file.exists(cif) && file.exists(rsif) && file.exists(dat)) ) {
          m_engine <- MutationEngine(setup_code = setup_ProCESS_code,
                                     COSMIC_account=COSMIC_account,
                                     context_sampling=c_sampling,
                                     max_repetition_storage=max_r_storage)
        }
      } else {
        file.remove(c(cif, rsif, dat))
      }

      measure <- peakRAM(
        m_engine <- MutationEngine(setup_code = setup_ProCESS_code,
                                   context_sampling=c_sampling,
                                   max_repetition_storage=max_r_storage))

      result <- data.frame(elapsed_time=measure$Elapsed_Time_sec,
                           load_saved = l_saved, repetition_id = r_id,
                           context_sampling = c_sampling,
                           max_repetition_storage = max_r_storage)
      
      warmup_counter <- warmup_counter-1
    }

    warmup_counter <- 1

    if (saving) {
      if (is.null(data_result)) {
        data_result <- result
      } else {
        data_result <- rbind(data_result, result)
      }
      saveRDS(data_result, output_file)
    }
  }

  current_test <- current_test + 1
}

plot1 <- ggplot(data_result %>% filter(context_sampling == 20),
                  aes(x = max_repetition_storage, y = elapsed_time/60,
                      color = load_saved, group = load_saved)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 3e+5) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Max. repetition storage",
    y = "Building time (minutes)"
  ) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  theme_minimal() +
  scale_color_discrete(
    name = "Index source", 
    labels = c("Built",
               "Loaded")
  )

save_image(section, "max_rep_storage-vs-time", extension=image_ext,
           directory=image_dir)

if (print_table) {
  cat("Max repetition storage-vs-Time\n")
  for (load_saved_v in c(TRUE, FALSE)) {
    if (load_saved_v == TRUE) {
      type <- "Loaded"
    } else {
      type <- "Built"
    }
    results <- compare_complexity(data_result %>%
                                    filter(context_sampling == 20 &
                                             load_saved == load_saved_v),
                                  "max_repetition_storage", "elapsed_time")
    latex_AICc_row(results$AICc, AICc_delta_threshold, first_columns=c(type))
  }
}

plot2 <- ggplot(data_result %>% filter(max_repetition_storage == 4e+5),
                aes(x = context_sampling, y = elapsed_time/60,
                    color = load_saved, group = load_saved)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 15) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  labs(
    x = "Context sampling",
    y = "Building time (minutes)"
  ) +
  theme_minimal() +
  scale_color_discrete(
    name = "Index source", 
    labels = c("Built",
               "Loaded")
  )

save_image(section, "context_sampling-vs-time", extension=image_ext,
           directory=image_dir)

if (print_table) {
  cat("\nContext sampling-vs-Time\n")
  for (load_saved_v in c(TRUE, FALSE)) {
    if (load_saved_v == TRUE) {
      type <- "Loaded"
    } else {
      type <- "Built"
    }
    results <- compare_complexity(data_result %>%
                                    filter(max_repetition_storage == 4e+5 &
                                             load_saved == load_saved_v),
                                  "context_sampling", "elapsed_time")
    latex_AICc_row(results$AICc, AICc_delta_threshold, first_columns=c(type))
  }
}


shared_legend <- get_legend(
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

y_limit <- range(data_result$elapsed_time/60)

plot1 <- plot1 + ylim(y_limit)
plot2 <- plot2 + ylim(y_limit)

plot1 <- plot1 + theme(legend.position = "none")
plot2 <- plot2 + theme(legend.position = "none",
                       axis.title.y = element_blank())

plots_row <- plot_grid(plot1, plot2, labels = c('A', 'B'), ncol = 2,
                       label_size = 12)

final_plot <- plot_grid(plots_row, shared_legend, align = 'h', axis = 'b',
                        rel_widths = c(3, 0.4))

save_image(section, "overall", width=10, extension=image_ext,
           directory=image_dir)
