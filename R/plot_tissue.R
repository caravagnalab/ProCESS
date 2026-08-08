## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2026 - Giulio Caravagna <gcaravagna@units.it>
##                           Alberto Casagrande <alberto.casagrande@uniud.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Plot a tissue
#'
#' @description This function plots the tissue
#' @details
#' This function represents cells distribution over a tissue. Each cells is
#' labelled and colored according to its label (see parameter
#' `label_function`). The tissue is draws as a heatmap of hexagonal bins for
#' efficiency porpoise.
#' @usage plot_tissue(simulation)
#' @param simulation A simulation object.
#' @param num_of_bins The number of bins (default: 100).
#' @param before_sample A sample name. When provided, this function represents
#'   the tissue as appeared just before the specified sampling. The parameters
#'   `before_sample` and `at_sample` are mutually exclusive (optional).
#' @param at_sample A sample name. When provided, this function represents
#'   the tissue as appeared when the specified sampling was about to be
#'   collected. The parameters `before_sample` and `at_sample` are mutually
#'   exclusive (optional).
#' @param plot_next_sample_regions A Boolean value. When `before_sample` is
#'   set and `plot_next_sample_regions` is set to be `TRUE`, this function
#'   plots the regions of the samples collected at the same simulated time
#'   of the specified sample. When, instead, `at_sample` is set and
#'   `plot_next_sample_regions` is set to be `TRUE`, the function plots the
#'   regions of the samples collected at the same simulated time of the
#'   specified sample, but not before the specified sample (default: `FALSE`).
#' @param plot_sample_region A Boolean value. When either `at_sample` or
#'   `before_sample` are set and `plot_sample_region` is set to be `TRUE`,
#'   the function also plots the region of the specified sample
#'   (default: `TRUE`).
#' @param color_map A named vector representing the color of the labels
#'   map (optional when `label_function` is `NULL`; mandatory, otherwise).
#' @param list_all_labels A Boolean flag to show all labels in
#'   the legend (default: `FALSE`).
#' @param label_function A function whose input is the result of
#'   <code>[TissueSimulation$get_cells()]</code> and that returns
#'   a string vector whose length is the number of rows in the input data
#'   frame. The strings are the labels of the corresponding cells and
#'   the function represents the different labels in the returned plot
#'   by coloring the cells according to `color_map`. If `label_function`
#'   is specified, then `color_map` becomes mandatory. When the parameter
#'   is set to `NULL`, the cells are labelled by their species names
#'   (default: `NULL`).
#' @param focus_function A function whose input is the result of
#'   <code>[TissueSimulation$get_cells()]</code> and that returns
#'   a Boolean vector whose length is the number of rows in the input data
#'   frame. When one the row in the output is `FALSE` the corresponding
#'   cells is plotted in grey. When the parameter is set to `NULL`, all
#'   tumour simulation cells are colored (default: `NULL`).
#' @param alpha_function A function whose input is the result of
#'   <code>[TissueSimulation$get_cells()]</code> and that returns
#'   a real vector whose values are in the interval \eqn{[0,1]} and whose
#'   length is the number of rows in the input data frame. Each value in
#'   the output is used as alpha level of the corresponding cell. When the
#'   parameter is set to `NULL`, all tumour simulation cells have alpha
#'   level `1` (default: `NULL`).
#' @return An editable ggplot plot.
#' @examples
#' # set the seed
#' set.seed(0)
#'
#' # build a tissue simulation
#' sim <- TissueSimulation(width = 600, height = 600)
#'
#' # avoid drift
#' sim$death_activation_level <- 50
#'
#' # add the mutant A
#' sim$add_mutant("A", c(duplication = 0.12, death = 0.05))
#'
#' # place a cell in the tissue and simulate it until 10 cells
#' sim$place_cell("A", 300, 300)
#' sim$run_up_to_size("A", 10, quiet = TRUE)
#'
#' # add the mutant B and let mutate a border cell of A in B
#' sim$add_mutant("B", c(duplication = 0.145, death = 0.06))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "B")
#'
#' # simulate the tissue up to 30 cells in B
#' sim$run_up_to_size("B", 30, quiet = TRUE)
#'
#' # add the third mutant and let one cell of A mutate into C
#' sim$add_mutant("C", c(duplication = 0.15, death = 0.06))
#' sim$mutate_progeny(sim$choose_border_cell_in("A"), "C")
#'
#' # simulate the tissue until C consists of 25000 cells
#' sim$run_up_to_size("C", 25000, quiet = TRUE)
#'
#' # collect the sample "S1"
#' sim$sample_cells("S1", c(145, 230), c(215, 300))
#'
#' # let the simulation reach 25000 cells in C again
#' sim$run_up_to_size("C", 25000, quiet = TRUE)
#'
#' # collect two samples
#' sim$sample_cells("S2", c(350, 300), c(420, 370))
#' sim$sample_cells("S3", c(200, 350), c(270, 420))
#'
#' # add a further mutant and derive it from B
#' sim$add_mutant(name = "D", c(duplication = 0.8, death = 0.01))
#' sim$mutate_progeny(sim$choose_border_cell_in("B"), "D")
#'
#' # let the tumour evolve until the mutant C and D cumulatively
#' # consist of 10000 cells
#' sim$run_until(sim$var("C") + sim$var("D") == 1e5, quiet = TRUE)
#'
#' # plot the tissue in the current status
#' plot_tissue(sim)
#'
#' # plot the tissue as it was when "S3" was about to be sampled
#' plot_tissue(sim, at_sample = "S3")
#'
#' # plot the tissue as it was when "S3" was about to be sampled and
#' # highlight the regions of the samples collected at the same
#' # simulated time, but not before it, i.e., "S3"
#' plot_tissue(sim, at_sample="S3",
#'             plot_next_sample_regions = TRUE)
#'
#' # plot the tissue as it was just before sampling "S3"
#' plot_tissue(sim, before_sample="S3")
#'
#' # plot the tissue as it was just before sampling "S3" and highlight
#' # the regions of the samples collected at the same simulated time,
#' # i.e., "S2" and "S3"
#' plot_tissue(sim, before_sample="S3",
#'             plot_next_sample_regions = TRUE)
#'
#' # define a custom color map
#' color_map <- c(A="#B2DF8A", B="#E31A1C", C="#C41E4E", D="#FEAAAA")
#' names(color_map) <- c("A", "B", "C", "D")
#'
#' plot_tissue(sim, color_map = color_map)
#'
#' # this function returns `TRUE` for cells in the rectangle
#' # [200,400]x[300,500]
#' focus_function <- function(cells) {
#'   (cells$position_x >= 200 & cells$position_x <= 400
#'    & cells$position_y >= 300 & cells$position_y <= 500)
#' }
#'
#' # plot the tissue highlighting the region [200,400]x[300,500]
#' plot_tissue(sim, focus_function = focus_function)
#'
#' # this function labels cells. Inside the rectangle
#' # [200,400]x[300,500] the label is the cell's mutant name.
#' # Outside, the rectangle label is cell's mutant name with
#' # "outside" attached.
#' label_function <- function(cells) {
#'   library(dplyr)
#'
#'   cells %>%
#'     mutate(label = if_else(cells$position_x >= 200
#'                            & cells$position_x <= 400
#'                            & cells$position_y >= 300
#'                            & cells$position_y <= 500,
#'                            .data$mutant,
#'                            paste(.data$mutant, "outside"))) %>%
#'     pull(label)
#' }
#'
#' # get the plot labels (i.e., mutants + paste(mutants, "outside"))
#' mutants <- sim$get_mutants() %>% dplyr::pull(mutant)
#' labels <- c(mutants, paste(mutants, "outside"))
#'
#' # create a color map for the labels
#' color_map <- RColorBrewer::brewer.pal(n = length(labels), name = "Set1")
#' names(color_map) <- labels
#'
#' # plot the tissue labelling the cells according to
#' # `label_function`. The parameter `color_map` is mandatory.
#' plot_tissue(sim, label_function = label_function, color_map = color_map)
#' @seealso [build_snapshot_video()]
#' @export
#'
plot_tissue <- function(simulation, num_of_bins = 100,
                        before_sample = NULL,
                        at_sample = NULL,
                        plot_next_sample_regions = FALSE,
                        plot_sample_region = TRUE,
                        color_map = NULL,
                        list_all_labels = FALSE,
                        label_function = NULL,
                        focus_function = NULL,
                        alpha_function = NULL) {
  stopifnot(inherits(simulation, "Rcpp_TissueSimulation"))

  sample_info <- NULL

  # Get the cells
  if (is.null(at_sample)) {
    if (is.null(before_sample)) {
      cells <- simulation$get_cells()
    } else {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% dplyr::filter(.data$name == before_sample)

      selected_samples <- samples %>%
        dplyr::filter(.data$time == sample_info$time[[1]])

      first_sample_at_time <- samples %>%
        dplyr::filter(.data$id == min(selected_samples$id))

      cells <- simulation$get_cells(first_sample_at_time$name[[1]])

      from_id <- first_sample_at_time$id[[1]]
    }
  } else {
    if (is.null(before_sample)) {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% dplyr::filter(.data$name == at_sample)

      from_id <- sample_info$id[[1]]

      cells <- simulation$get_cells(at_sample)
    } else {
      stop(paste("The parameters `at_sample` and `before_sample`",
                 "are mutually exclusive."))
    }
  }

  if (is.null(label_function)) {
    cells <- cells %>% add_species_col("label")

    if (is.null(color_map)) {
      color_map <- get_species_colors(simulation)
    }
  } else {
    if (is.null(color_map)) {
      stop("\"color_map\" is mandatory when \"label_function\" ",
           "is specified.")
    }

    cells[["label"]] <- label_function(cells)
  }

  if (is.null(alpha_function)) {
    cells$alpha_level <- 1
  } else {
    cells$alpha_level <- alpha_function(cells)
  }

  cells[["label"]] <- factor(cells[["label"]],
                             levels = unique(names(color_map)))

  if (is.null(focus_function)) {
    highlighted_cells <- cells
    masked_cells <- cells[0, ]
  } else {
    highlighted_cells <- cells %>%
      dplyr::filter(focus_function(cells))
    masked_cells <- cells %>%
      dplyr::filter(!focus_function(cells)) %>%
      dplyr::mutate(alpha_level = .data$alpha_level * 0.7)
  }

  tissue_sizes <- simulation$get_tissue_size()

  bin_width_x <- tissue_sizes[1] / num_of_bins
  bin_width_y <- tissue_sizes[2] / num_of_bins

  pl <- ggplot2::ggplot() +
    ggplot2::stat_summary_hex(
      data = masked_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   fill = ggplot2::stage(start = .data$label,
                                         after_stat = NULL),
                   z = .data$alpha_level,
                   alpha = ggplot2::after_stat(value)),
      fun = max,
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = FALSE
    ) +
    ggplot2::stat_summary_hex(
      data = masked_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   z = .data$alpha_level,
                   fill = NA,
                   alpha = ggplot2::after_stat(value)),
      fill = "grey70",
      fun = max,
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = FALSE
    ) +
    ggplot2::stat_summary_hex(
      data = highlighted_cells,
      ggplot2::aes(x = .data$position_x,
                   y = .data$position_y,
                   fill = ggplot2::stage(start = .data$label,
                                         after_stat = NULL),
                   z = .data$alpha_level,
                   alpha = ggplot2::after_stat(value)),
      fun = max,
      binwidth = c(bin_width_x, bin_width_y),
      show.legend = TRUE
    ) +
    ggplot2::guides(alpha = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_fill_manual(values = color_map,
                               limits = names(color_map),
                               drop = !list_all_labels) +
    my_theme() +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlim(-10, tissue_sizes[1] + 10) +
    ggplot2::ylim(-10, tissue_sizes[2] + 10)

  if (!is.null(sample_info)) {
    if (plot_sample_region) {
      pl <- pl + ggplot2::annotate("rect",
                                   xmin = sample_info$xmin[[1]],
                                   xmax = sample_info$xmax[[1]],
                                   ymin = sample_info$ymin[[1]],
                                   ymax = sample_info$ymax[[1]],
                                   fill = NA, color = "black")
    }

    if (plot_next_sample_regions) {
      if (is.null(before_sample)) {
        selected_samples <- samples %>%
          dplyr::filter(.data$time == sample_info$time[[1]],
                        .data$id >= from_id)
      }

      for (row_idx in seq_len(nrow(selected_samples))) {
        sample_info <- selected_samples[row_idx, ]
        pl <- pl + ggplot2::annotate("rect",
                                     xmin = sample_info$xmin[[1]],
                                     xmax = sample_info$xmax[[1]],
                                     ymin = sample_info$ymin[[1]],
                                     ymax = sample_info$ymax[[1]],
                                     fill = NA, color = "black")
      }
    }
  }

  pl
}

#' Building a video of the snapshots
#'
#' @description This function builds a video of the snapshots.
#' @details This function builds a video of the snapshots. It is available
#'   only if the package `av` is installed.
#' @usage build_snapshot_video(simulation)
#' @param simulation A simulation object.
#' @param output_file The path of the output video. When it is set to
#'   `NULL`, the output video path has the format
#'   `<simulation path>_evolution.mp4` (default: `NULL`).
#' @param plot_function The function used to plot each frame of the video.
#'   When is set to `NULL`, the function `plot_tissue()` is used (default:
#'   `NULL`).
#' @param width The width of the video (default: `800`).
#' @param height The height of the video (default: `600`).
#' @param framerate The video framerate in frame/sec (default: `1`).
#' @param res The video resolution (default: `150`).
#' @param quiet A Boolean flag to enable/disable the messages (default:
#'   `FALSE`).
#' @param workers The number of parallel processes generating frames.
#'   This parameter is used only when the packages "furrr" and
#'   "progressr" are installed. When it is set to `NULL`, the function
#'   uses as many processes as the number of available processors minus
#'   one (default: `NULL`).
#' @returns The name of the produced video file path.
#' @seealso `vignette("videos")`, [plot_tissue()]
#' @export
#'
build_snapshot_video <- function(simulation, output_file = NULL,
                                 plot_function = NULL,
                                 width = 800, height = 600,
                                 framerate = 1, res = 150,
                                 quiet = FALSE,
                                 workers = NULL) {

  if (!requireNamespace("av", quietly = TRUE)) {
    stop(
      "Package 'av' is required for using this function. ",
      "Please install it using install.packages('av').",
      call. = FALSE
    )
  }

  # if plot_function is NULL, use the plot_tissue function
  if (is.null(plot_function)) {
    # choose a color map that contains all the known species
    color_map <- get_species_colors(simulation)

    plot_function <- function(snapshot) {
      plot_tissue(snapshot,
                  color_map = color_map,
                  plot_sample_region = FALSE,
                  list_all_labels = TRUE)
    }
  }

  snapshot_files <- simulation$get_snapshot_info()[["file"]]

  if (is.null(output_file)) {
    output_file <- paste0(simulation$get_name(), "_evolution.mp4")
  }

  frame_dir <- tempfile(pattern = "ProCESS_video_")
  dir.create(frame_dir)

  on.exit(unlink(frame_dir, recursive = TRUE), add = TRUE)

  if (!quiet) {
    cat("Collecting frames...\n")
  }

  get_frame_filename <- function(order) {
    file.path(frame_dir,
              paste0("video_", sprintf("%010d", order), ".png"))
  }

  frame_files <- lapply(seq_along(snapshot_files), get_frame_filename)

  files_enum <- Map(function(ff, sf) {
    list(frame_file = ff, snapshot_file = sf)
  }, frame_files, snapshot_files)

  get_frames_par <- function(files_enum, workers, quiet) {
    # backup previous handlers and reset them on exit
    old_handlers <- progressr::handlers()
    on.exit(progressr::handlers(old_handlers), add = TRUE)

    if (requireNamespace("progress", quietly = TRUE)) {
      progressr::handlers(
        progressr::handler_progress(
          format = "[:bar] :current/:total (:percent) | ETA: :eta",
          clear = FALSE
        )
      )
    } else {
      progressr::handlers(
        progressr::handler_txtprogressbar(
          char = "=",
          style = 3
        )
      )
    }

    future::plan(future::multisession,
                 workers = min(workers, length(files_enum)))

    progressr::with_progress({
      pb <- progressr::progressor(along = files_enum)

      furrr::future_map(
        .x = files_enum,
        .f = function(files_enum_value) {
          snapshot <- recover_simulation(files_enum_value$snapshot_file,
                                         quiet = TRUE)

          plt <- plot_function(snapshot)

          ggplot2::ggsave(files_enum_value$frame_file, plot = plt,
                          width = width, height = height, dpi = res,
                          units = "px")

          # advance the progress bar
          pb()
        },
        .options = furrr::furrr_options(packages = c("dplyr", "ProCESS"),
                                        seed = TRUE, scheduling = FALSE)
      )
    }, enable = !quiet)

    future::plan(future::sequential)
  }

  get_frames_seq <- function(files_enum, quiet) {
    if (!quiet) {
      progress_bar <- txtProgressBar(min = 0, max = length(files_enum),
                                     style = 3)
    }
    idx <- 0
    if (!quiet) {
      setTxtProgressBar(progress_bar, idx)
    }
    for (files_enum_value in files_enum) {
      snapshot <- recover_simulation(files_enum_value$snapshot_file)

      plt <- plot_function(snapshot)

      ggplot2::ggsave(files_enum_value$frame_file, plot = plt,
                      width = width, height = height, dpi = res,
                      units = "px")
      if (!quiet) {
        idx <- idx + 1
        setTxtProgressBar(progress_bar, idx)
      }
    }

    if (!quiet) {
      close(progress_bar)
    }
  }

  if (requireNamespace("furrr", quietly = TRUE)
      && requireNamespace("progressr", quietly = TRUE)) {

    if (is.null(workers)) {
      workers <- future::availableCores() - 1
    }
    get_frames_par(files_enum, workers = workers, quiet)
  } else {
    if (!is.null(workers)) {
      warning(paste("Packages \"furrr\" and \"progressr\" are not available.",
                    "Parameter \"workers\" is not used."))
    }

    get_frames_seq(files_enum, quiet)
  }

  av::av_encode_video(
    input = unlist(frame_files),
    output = output_file,
    framerate = framerate,
    verbose = !quiet
  )
}