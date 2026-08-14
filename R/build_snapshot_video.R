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

normalize_pause <- function(x, framerate) {
  if (!"length" %in% names(x)) {
    stop("The parameter \"pauses_on_frame\" must be list of named list ",
         "`list(frame, length)`.")
  }

  if (inherits(x$length, "difftime")) {
    x$length <- as.numeric(x$length, units = "secs") * framerate

    return(x)
  }

  if (is.numeric(x$length) || is.integer(x$length)) {
    if (x$length != trunc(x$length) || x$length <= 0) {
      stop("When a \"pauses_on_frame\"' \"length\" value represent ",
           "the length's number of frames, it must be a ",
           "positive integer.")
    }

    return(x)
  }
  stop("The \"pauses_on_frame\"' \"length\" values must either be ",
       "the length's number of frames or a `difftime`.")
}

add_pauses <- function(frame_files, pauses, framerate) {

  if (is.null(pauses)) {
    pauses <- list()
  }

  frames <- tryCatch({
    sapply(pauses, function(x) x$frame)
  }, error = function(e) {
    stop("The \"pauses_on_frame\" parameter must be a list of named ",
         "lists `list(frame, length)`")
  })

  pauses <- pauses[order(frames)]

  pauses <- lapply(pauses, function(x) {
    normalize_pause(x, framerate = framerate)
  })

  new_frame_files <- list()

  pauses_idx <- 1
  for (frame_idx in seq_along(frame_files)) {
    new_frame_files <- append(new_frame_files,
                              frame_files[[frame_idx]])
    if (pauses_idx <= length(pauses)) {
      if (pauses[[pauses_idx]]$frame == frame_idx) {
        for (i in 1:pauses[[pauses_idx]]$length) {
          new_frame_files <- append(new_frame_files,
                                    frame_files[[frame_idx]])
        }
        pauses_idx <- pauses_idx + 1
      }
    }
  }

  new_frame_files
}

get_frames <- function(pause_length, framerate) {
  if (inherits(pause_length, "difftime")) {
    units(pause_length) <- "secs"
    pause_length <- pause_length * framerate
  }

  pause_length
}

get_event_pauses <- function(simulation, pauses_on_event, framerate) {

  pauses <- list()
  if (!is.null(pauses_on_event)) {
    if ("sampling" %in% names(pauses_on_event)) {
      if (inherits(pauses_on_event[["sampling"]], "difftime")
          || is.numeric(pauses_on_event[["sampling"]])) {
        samples <- simulation$get_samples_info()
        samples$length <- pauses_on_event[["sampling"]]
        sampling_pauses <- as.list(samples[["length"]])
        names(sampling_pauses) <- as.list(samples[["name"]])

        pauses_on_event[["sampling"]] <- sampling_pauses
      }
    }

    if ("mutant emerged" %in% names(pauses_on_event)) {
      if (inherits(pauses_on_event[["mutant emerged"]], "difftime")
          || is.numeric(pauses_on_event[["mutant emerged"]])) {
        mutants <- simulation$get_mutants()
        mutants$length <- pauses_on_event[["mutant emerged"]]
        mutant_pauses <- as.list(mutants[["length"]])
        names(mutant_pauses) <- mutants[["mutant"]]

        pauses_on_event[["mutant emerged"]] <- mutant_pauses
      }
    }

    snapshot_files <- simulation$get_snapshot_info()[["file"]]

    for (snapshot_idx in seq_along(snapshot_files)) {
      snapshot <- recover_simulation(snapshot_files[[snapshot_idx]],
                                     quiet = TRUE)

      new_events <- snapshot$get_just_occurred_events()

      pause_length <- 0
      for (event in list("sampling", "mutant emerged")) {
        if (event %in% names(pauses_on_event)) {
          if (event %in% names(new_events)) {
            event_pauses <- pauses_on_event[[event]]
            event_obj_name <- new_events[[event]]
            if (event_obj_name %in% names(event_pauses)) {
              pauses_frames <- get_frames(event_pauses[[event_obj_name]],
                                          framerate)
              pause_length <- max(pause_length, pauses_frames)
            }
          }
        }
      }

      if ("rate update" %in% names(pauses_on_event)) {
        if ("rate update" %in% new_events) {
          pauses_frames <- get_frames(pauses_on_event[["rate update"]],
                                      framerate)
          pause_length <- max(pause_length, pauses_frames)
        }
      }

      if (pause_length > 0) {
        pauses[[length(pauses) + 1]] <- list(frame = snapshot_idx,
                                             length = pause_length)
      }
    }
  }

  pauses
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
#' @param pauses_on_event A named list specifying the pauses on event in
#'   the output video.
#'   The names of the list must be among "mutant rising" and "sampling".
#'   The values can either be a numeric value, a `difftime` object, or
#'   a named list. The numeric value represents the pause length in number
#'   of frames. The `difftime` object denotes the pause length. Finally,
#'   the named list described the pauses for specific events. When the name
#'   of the element is "mutant rising", the names of the sub-list are among
#'   the simulated mutants. When, instead, the name of the element is
#'   "sampling", the names of the sub-list are among the collected samples.
#'   In both the cases, the values represent the pause lengths for the
#'   specific event and can be either a numeric value or a `difftime`
#'   object as for the generic specification. When `pauses_on_event` is
#'   `NULL`, no pauses are added to the output video (default: `NULL`).
#' @param pauses_on_frame A list specifying the pauses on frame in the
#'   output video.
#'   Each element of the list is a named list whose names are `frame` and
#'   `length`. The element `frame` is the video frame from which the pause
#'   is aimed. Instead, `length` is the pause's length expressed in either
#'   number of frames, when `length` is a numeric value, or a `difftime`.
#'   When `pauses_on_frame` is NULL, no pauses are added to the output video
#'   (default: `NULL`).
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
                                 pauses_on_event = NULL,
                                 pauses_on_frame = NULL,
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

  if (!is.null(pauses_on_event)) {
    pauses_on <- get_event_pauses(simulation, pauses_on_event, framerate)
    frame_files <- add_pauses(frame_files, pauses_on, framerate)
  }

  if (!is.null(pauses_on_frame)) {
    frame_files <- add_pauses(frame_files, pauses_on_frame, framerate)
  }

  av::av_encode_video(
    input = unlist(frame_files),
    output = output_file,
    framerate = framerate,
    verbose = !quiet
  )
}