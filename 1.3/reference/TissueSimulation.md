# Simulating the cell evolution in a tissue

This class simulates the cell evolution on a tissue.

This method builds a new simulation.

## Usage

``` r
TissueSimulation(name, width = 1000, height = 1000, save_snapshots = FALSE)
```

## Arguments

- name:

  The name of the simulation (default:
  "`ProCESS_<year>_<hour><minute><second>`").

- width:

  The width of the simulated tissue (default: 1000).

- height:

  The height of the simulated tissue (default: 1000).

- save_snapshots:

  A flag to save simulation snapshots on disk (default: `FALSE`).

- rates:

  A data frame specifying the simulation species and their rates
  (default: `NULL`). See
  [`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md)
  for the data frame specification. Differently from
  [`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md)
  the species are automatically added to the new simulation.

- epigenetic_states:

  A list of epigenetic states (default: `NULL`).

- seed:

  The seed for the pseudo-random generator (optional).

## Details

The objects of this class can simulate the evolution of many cells
belonging to different *species* on a tissue. Each cell can duplicate or
die according to the rates that delineate the cell species.

`TissueSimulation()` supports epigenetic evolutions, and it lets users
define species pairs that belong to the same mutant (even though, its
genomic characterisation is unknown) and differ because of their
epigenetic states.

`TissueSimulation()` models epigenetic mutations and allows a cell in
one of mutant species to generate a new cell belonging to the other
species of the same mutant at a specified rate.

`TissueSimulation()` also allows users to schedule mutations from one
mutant to a different mutant.

## Examples

``` r
# create a TissueSimulation object storing binary dump in a temporary
# directory. The data are deleted from the disk as soon as the object
# is destroyed.
sim <- TissueSimulation("test")

# the name of the simulation is "test"
sim$get_name()
#> [1] "test"

# however no directory "test" has been created in the working directory
"test" %in% list.files()
#> [1] FALSE

# By using the optional parameter `save_snapshots`, we force the
# simulation to save its progresses in a local directory whose name
# is the name of the simulation, i.e., "test". This data will be
# preserved when the simulation object will be destroyed.
sim <- TissueSimulation("test", save_snapshots = TRUE)

# the directory "test" exists and contains a binary dump of
# the simulation
"test" %in% list.files()
#> [1] TRUE

# the directory persists even after the object destruction
rm(sim)
"test" %in% list.files()
#> [1] TRUE

# let us manually delete the "test" directory
unlink("test", recursive = TRUE)

# the name parameter is optional
sim <- TissueSimulation(save_snapshots = TRUE)

# the name of the simulation is `ProCESS_<YY><MM><DD>_<HH><MM><SS>`
sim$get_name()
#> [1] "ProCESS_20260709-003928"

# the simulation dump have been saved in a directory named
# after the simulation name
list.files(pattern = "^ProCESS_")
#> [1] "ProCESS_20260709-003928"

# let us remove the object and manually delete the simulation
# directory
rm(sim)
unlink(list.files(pattern = "^ProCESS_"), recursive = TRUE)

# users can provide a random seed to the simulation...
sim <- TissueSimulation(seed = 13)

# ..., specify the size of the simulated space by using the
# optional parameters `width` and `height`, or...
sim <- TissueSimulation(width = 1200, height = 900)
sim$get_tissue_size()
#> [1] 1200  900

# ... build a simulation, add its species, and set their rates
# by passing a data frame specifying the event rates
df_rates <- data.frame(
  mutant = c("A", "A", "B"),
  event = c("duplication", "death", "duplication"),
  rate = c(0.3, 0.2, 0.2)
)
df_rates
#>   mutant       event rate
#> 1      A duplication  0.3
#> 2      A       death  0.2
#> 3      B duplication  0.2

sim <- TissueSimulation(rates = df_rates)
sim
#> ──  ProCESS   D   S   M  ProCESS_20260709-003928 ────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ── Species: 2, without epigenetics 
#>    
#>    =======  ===  ===  ======  ===
#>    species   λ    δ   counts   % 
#>    =======  ===  ===  ======  ===
#>          A  0.3  0.2    0     NaN
#>          B  0.2  0.0    0     NaN
#>    =======  ===  ===  ======  ===
#> 
#> ── Firings: 0 total 
#> 
#> ✖ The simulation has no samples yet!

# if epigenetic states are needed, the data frame must also contain
# the columns `epistate` and `first.child.epistate`
df_rates <- data.frame(
  mutant = c("A", "A", "A", "A", "A", "B", "B"),
  epistate = c("E1", "E1", "E1", "E1", "E2", "E1", "E1"),
  event = c("duplication", "death", "switch", "switch",
            "duplication", "duplication", "switch"),
  first.child.epistate = c("E1", NA, "E2", "E3", "E2", "E1", "E3"),
  rate = c(0.3, 0.2, 0.01, 0.04, 0.2, 0.2, 0.01)
)
df_rates
#>   mutant epistate       event first.child.epistate rate
#> 1      A       E1 duplication                   E1 0.30
#> 2      A       E1       death                 <NA> 0.20
#> 3      A       E1      switch                   E2 0.01
#> 4      A       E1      switch                   E3 0.04
#> 5      A       E2 duplication                   E2 0.20
#> 6      B       E1 duplication                   E1 0.20
#> 7      B       E1      switch                   E3 0.01

sim <- TissueSimulation(rates = df_rates)
sim
#> ──  ProCESS   D   S   M  ProCESS_20260709-003928 ────────────────────────────────────── ▣  [1000x1000]  ⏱ 0 ──
#> 
#> ── Species: 6, with epigenetics 
#>    
#>    =======  ===  ===  ======  ===
#>    species   λ    δ   counts   % 
#>    =======  ===  ===  ======  ===
#>      A[E1]  0.3  0.2    0     NaN
#>      A[E2]  0.2  0.0    0     NaN
#>      A[E3]  0.0  0.0    0     NaN
#>      B[E1]  0.2  0.0    0     NaN
#>      B[E2]  0.0  0.0    0     NaN
#>      B[E3]  0.0  0.0    0     NaN
#>    =======  ===  ===  ======  ===
#> 
#> ── Epigenetic switches 
#>    
#>    =======  ====  =====
#>    species   ε     dest
#>    =======  ====  =====
#>      A[E1]  0.01  A[E2]
#>      A[E1]  0.04  A[E3]
#>      B[E1]  0.01  B[E3]
#>    =======  ====  =====
#> 
#> ── Firings: 0 total 
#> 
#> ✖ The simulation has no samples yet!
```
