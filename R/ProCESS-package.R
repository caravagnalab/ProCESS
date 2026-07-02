## This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
## Copyright (C) 2023-2025 - Alberto Casagrande <alberto.casagrande@uniud.it>
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


.onLoad <- function(libname, pkgname) {
  loadModule("Mutants", TRUE)
  loadModule("Mutations", TRUE)
  loadModule("Sequencing", TRUE)
  loadModule("Logics", TRUE)

  if (exists("WGD", envir = .GlobalEnv, inherits = FALSE)) {
    message(paste(pkgname, "overwrites active binding \"WGD\"."))

    rm("WGD", envir = .GlobalEnv)
  }

  wg_doubling <- function() new(WholeGenomeDoubling)

  makeActiveBinding("WGD", wg_doubling, .GlobalEnv)

  ## usethis namespace: start
  ## usethis namespace: end
    NULL

}

.onAttach <- function(libname, pkgname) {
  setMethod("+", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)
  setMethod("+", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_sum(e1, e2)
            }, where = .GlobalEnv)

  setMethod("-", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)
  setMethod("-", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_subtract(e1, e2)
            }, where = .GlobalEnv)

  setMethod("*", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)
  setMethod("*", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_multiply(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "numeric", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Expression", e2 = "numeric"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "numeric", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Variable", e2 = "numeric"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Expression", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Expression"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod(">", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_gt(e1, e2)
            }, where = .GlobalEnv)
  setMethod(">=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ge(e1, e2)
            }, where = .GlobalEnv)
  setMethod("==", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_eq(e1, e2)
            }, where = .GlobalEnv)
  setMethod("!=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_ne(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<=", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_le(e1, e2)
            }, where = .GlobalEnv)
  setMethod("<", signature(e1 = "Rcpp_Variable", e2 = "Rcpp_Variable"),
            function(e1, e2) {
              logics_lt(e1, e2)
            }, where = .GlobalEnv)

  setMethod("&", signature(e1 = "Rcpp_Formula", e2 = "Rcpp_Formula"),
            function(e1, e2) {
              logics_and(e1, e2)
            }, where = .GlobalEnv)
  setMethod("|", signature(e1 = "Rcpp_Formula", e2 = "Rcpp_Formula"),
            function(e1, e2) {
              logics_or(e1, e2)
            }, where = .GlobalEnv)

  #' @importFrom rlang .data
  setMethod("show", "Rcpp_TissueSimulation", function(object) {
    # If it can be simulated
    sim_status <- nrow(object$get_cells())

    # If it has samples assigned
    sam_status <- sim_status & nrow(object$get_samples_info())

    mut_status <- FALSE

    sim_status <- ifelse(
      sim_status,
      crayon::bgGreen(crayon::white(" D ")),
      crayon::bgRed(crayon::white(" D "))
    )

    sam_status <- ifelse(
      sam_status,
      crayon::bgGreen(crayon::white(" S ")),
      crayon::bgRed(crayon::white(" S "))
    )

    mut_status <- ifelse(
      mut_status,
      crayon::bgGreen(crayon::white(" M ")),
      crayon::bgRed(crayon::white(" M "))
    )

    cli::cli_rule(
      left = paste(
        crayon::bgYellow(" ProCESS "),
        sim_status, sam_status, mut_status,
        object$get_name()
      ),
      right = paste0(
        crayon::red("\u25A3 "),
        " [",
        paste(object$get_tissue_size(), collapse = "x"),
        "]",
        crayon::red("  \u23f1 "),
        signif(object$get_clock(), digits = 3)
      )
    )

    species_name <- function(mutant, epistate) {
      dplyr::if_else(is.na(epistate),
                     epistate,
                     paste0(mutant, "[", epistate, "]"))
    }

    counts_tab <- object$get_counts()
    if (nrow(counts_tab) > 0) {
      counts_tab <- counts_tab %>%
        dplyr::mutate(
          `%` = 100 * .data$counts / sum(.data$counts)
        )

      rates <- object$get_rates()

      with_epigenetics <- "epistate" %in% colnames(rates)
      if (with_epigenetics) {
        with_epigenetics <- max(nchar(rates$epistate)) != 0
      }

      if (with_epigenetics) {
        rates <- rates %>%
          dplyr::mutate(species = species_name(.data$mutant, .data$epistate),
                        dest = species_name(.data$mutant,
                                            .data$first.child.epistate)) %>%
          dplyr::select(.data$species, .data$event,
                        .data$dest, .data$rate)
        counts_tab <- counts_tab %>%
          dplyr::mutate(species = species_name(.data$mutant,
                                               .data$epistate)) %>%
          dplyr::select(.data$species, .data$counts, .data$`%`)
      } else {
        rates <- rates %>%
          dplyr::mutate(species = .data$mutant) %>%
          dplyr::select(.data$species, .data$event, .data$rate)
        counts_tab <- counts_tab %>%
          dplyr::mutate(species = .data$mutant) %>%
          dplyr::select(.data$species, .data$counts, .data$`%`)
      }

      death_rates <- rates %>% dplyr::filter(.data$event == "death") %>%
        dplyr::mutate(death = .data$rate) %>%
        dplyr::select(.data$species, .data$death)

      duplication_rates <- rates %>%
        dplyr::filter(.data$event == "duplication") %>%
        dplyr::mutate(duplication = .data$rate) %>%
        dplyr::select(.data$species, .data$duplication)

      my_tab <- duplication_rates %>%
        dplyr::full_join(death_rates, by = c("species")) %>%
        dplyr::full_join(counts_tab, by = c("species")) %>%
        dplyr::mutate(death = dplyr::coalesce(.data$death, 0),
                      duplication = dplyr::coalesce(.data$duplication, 0)) %>%
        dplyr::select(.data$species, .data$duplication,
                      .data$death, .data$counts, .data$`%`) %>%
        dplyr::arrange(.data$species)
      rm(duplication_rates)
      rm(death_rates)

      colnames(my_tab)[2:3] <- c(" \u03BB ", " \u03B4 ") #, " \u03B5 ")

      cli::cli_h3(text = paste0("Species: {.field {nrow(my_tab)}}, ",
                                "{.field {ifelse(with_epigenetics, ",
                                "crayon::green('with'), ",
                                "crayon::red('without'))}} epigenetics"))

      cat("   ", knitr::kable(my_tab, format = "rst", align = "rcccc"),
          sep = "\n   ")

      if (with_epigenetics) {

        cli::cli_h3(text = "Epigenetic switches")
        switch_rates <- rates %>%
          dplyr::filter(.data$event != "death" & !is.na(.data$dest) &
                        .data$species != .data$dest) %>%
          dplyr::select(.data$species, .data$rate, .data$dest) %>%
          dplyr::arrange(.data$species, .data$dest)

        colnames(switch_rates)[2] <- " \u03B5 "

        cat("   ", knitr::kable(switch_rates, format = "rst", align = "rcr"),
            sep = "\n   ")

        rm(switch_rates)
      }

      rm(rates)

      f_table <- object$get_firings()

      if (with_epigenetics) {
        f_table <- f_table %>%
          dplyr::mutate(species = species_name(.data$mutant, .data$epistate))
      } else {
        f_table <- f_table %>% dplyr::mutate(species = .data$mutant)
      }

      cli::cli_h3(text = paste("Firings:", sum(f_table$fired), "total"))

      if (nrow(f_table) > 0) {
        nch <- f_table$fired %>% nchar %>% max

        for (s in my_tab$species) {
          s_ftab <- f_table %>% dplyr::filter(species == s) %>%
            dplyr::arrange(.data$event)

          if (with_epigenetics) {
            sprintf(
              paste0("\n\tSpecies %s: %", nch, "s (deaths), %",
                     nch, "s (duplications) and %", nch, "s (switches)"),
              s,
              s_ftab$fired[1],
              s_ftab$fired[2],
              s_ftab$fired[3]
            ) %>% cat()
          } else {
            sprintf(
              paste0("\n\tSpecies %s: %", nch, "s (deaths) and %",
                     nch, "s (duplications)"),
              s,
              s_ftab$fired[1],
              s_ftab$fired[2]
            ) %>% cat()
          }
        }
      }
    }


    samples <- object$get_samples_info()
    if (nrow(samples) == 0) {
      cat("\n")
      cli::cli_alert_danger("The simulation has no samples yet!")
    } else {
      t_samples <- samples$time %>% unique()

      is_longitudinal <- length(t_samples) > 1
      is_multiregion <- any(table(samples$time) > 1)

      is_longitudinal <- ifelse(is_longitudinal,
                                crayon::bgGreen(crayon::white(" Yes ")),
                                crayon::bgRed(crayon::black(" No ")))
      is_multiregion <- ifelse(is_multiregion,
                               crayon::bgGreen(crayon::white(" Yes ")),
                               crayon::bgRed(crayon::black(" No ")))

      cat("\n")

      ncells <- samples$tumour_cells %>% sum

      cli::cli_h3(text = paste0("Samples ({.emph multi-region} ",
                                "{is_multiregion}, {.emph longitudinal} ",
                                "{is_longitudinal}): {.field {ncells}} ",
                                "cells from {.field {nrow(samples)}} ",
                                "samples at {.field {length(t_samples)}} ",
                                "timepoints"))

      cat("   ",
          knitr::kable(samples, format = "rst", align = "rcrccc"),
          sep = "\n   ")

    }
  }, where = .GlobalEnv)
}
