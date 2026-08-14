# Scheduling a mutation

This method schedules a mutant mutation

## Arguments

- src:

  The name of the mutant from which the mutation occurs.

- dest:

  The name of the mutant to which the mutation leads.

- time:

  The simulated time at which the mutation will occurs.

## Details

The mutation can occur from any of the species of the source mutant to
the species of the destination mutant with a consistent epigenetic
state. For the sake of example, if the simulation has no epigenetic
states and a mutation from "A" to "B" is scheduled, then the first
duplication of an "A"'s cell occurring after the specified time
generates two cells: one of them belong to "A" and the other to "B".
Analogously, if the simulation has epigenetic states and the first
duplication of an "A"'s cell after the specified time occurs to a cell
in the species "AEi", then the offspring of the cell consists of one
cell in "AEi" and one in "BEi".

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# add mutant "B" and set its rates
sim$add_mutant("B", c(duplication = 0.4, death = 0.01))

# schedule an evolution from mutant "A" to mutant "B" at time 50
sim$schedule_mutation(src = "A", dst = "B", time = 50)

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation up to the first cell in B
sim$run_up_to_size("B", 1)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


sim
#> ──  ProCESS   D   S   M  ProCESS_20260814-004521 ─────── ▣  [1000x1000]  ⏱ 50 ──
#> 
#> ── Species: 2, without epigenetics 
#>    
#>    =======  ===  ====  ======  ==========
#>    species   λ    δ    counts      %     
#>    =======  ===  ====  ======  ==========
#>          A  0.2  0.01   593    99.8316498
#>          B  0.4  0.01    1     0.1683502 
#>    =======  ===  ====  ======  ==========
#> 
#> ── Firings: 731 total 
#> 
#>  Species A:  69 (deaths) and 662 (duplications)
#>  Species B:   0 (deaths) and   0 (duplications)
#> ✖ The simulation has no samples yet!
```
