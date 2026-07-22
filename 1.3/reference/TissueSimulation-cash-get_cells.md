# Getting the tissue cells

This method returns information about tumour tissue cells

## Arguments

- lower_corner:

  The lower-left corner of the selection frame (optional).

- upper_corner:

  The upper-right corner of the selection frame (optional).

- mutant_filter:

  The vector of the to-be-selected mutant names (optional).

- epigenetic_filter:

  The vector of the to-be-selected epigenetic states (optional).

## Value

A data frame having the columns `cell_id`, `mutant`, `position_x`,
`position_y`, and, when the simulation has epigenetic states,
`epistate`. Each row in the data frame describes a cell that satisfies
the provided filters and lays in the specified frame.

## Details

It collects some data about the cells in the tissue without altering the
tissue itself. The pairs of optional parameters `lower_corner` and
`upper_corner` define a frame of the tissue in which the data are
sampled. The optional parameters `mutant_filter` and `epigenetic_filter`
filter the collected cell data according to the cell mutant and
epigenetic state.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

sim$death_activation_level <- 100

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# add mutant "B" and set its species rates
sim$add_mutant("B",
               list(E1 = list(duplication = 0.3, death = 0.1, E2 = 0.02),
                    E2 = list(duplication = 0.1, death = 0.01, E1 = 0.01)))

# schedule a mutation from "A" to "B"
sim$schedule_mutation("A", "B", 20)

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                       


# collect all cells in the tissue
cells <- sim$get_cells()

# print some of them
head(cells)
#>   cell_id mutant epistate position_x position_y
#> 1    8079      B       E1        473        499
#> 2    8031      B       E1        473        500
#> 3    7959      B       E1        473        504
#> 4    7067      B       E1        473        505
#> 5    7822      B       E1        473        506
#> 6    8078      B       E1        474        499

# let us define a function to print some statistics
print_statistics <- function(cells) {
  cat(paste("num of cells:", nrow(cells)))
  cat(paste("\nmutants:", paste(unique(cells$mutant), collapse = " ")))
  cat(paste("\nepigenetic states:", paste(unique(cells$epistate),
                                          collapse = " ")))
  cat(paste0("\nframe: [",
             paste(range(cells$position_x), collapse = ","), "]x[",
             paste(range(cells$position_y), collapse = ","), "]\n\n"))
}

print_statistics(cells)
#> num of cells: 1735
#> mutants: B A
#> epigenetic states: E1 E2
#> frame: [473,522]x[479,532]
#> 

# get the cells in the frame [495,505]x[490,500]
cells <- sim$get_cells(lower_corner = c(495, 495),
                       upper_corner = c(505, 505))

print_statistics(cells)
#> num of cells: 118
#> mutants: A B
#> epigenetic states: E1 E2
#> frame: [495,505]x[495,505]
#> 

# cells can be filtered by mutant name...
cells <- sim$get_cells(mutant_filter = c("A"),
                       epigenetic_filter = c("E1", "E2"))

print_statistics(cells)
#> num of cells: 833
#> mutants: A
#> epigenetic states: E1 E2
#> frame: [475,522]x[479,523]
#> 

# ..., by epigenetic state, ...
cells <- sim$get_cells(mutant_filter = c("A", "B"),
                       epigenetic_filter = c("E2"))

print_statistics(cells)
#> num of cells: 460
#> mutants: B A
#> epigenetic states: E2
#> frame: [475,521]x[479,531]
#> 

# ..., or by position, mutant, and epigenetic state.
cells <- sim$get_cells(lower_corner = c(495, 495), upper_corner = c(505, 505),
                       mutant_filter = c("A", "B"),
                       epigenetic_filter = c("E2"))

print_statistics(cells)
#> num of cells: 45
#> mutants: B A
#> epigenetic states: E2
#> frame: [495,505]x[495,505]
#> 
```
