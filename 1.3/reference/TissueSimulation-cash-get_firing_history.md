# Getting the fired event history

This method returns a data frame reporting the number of events fired up
to each sampled simulation time.

## Value

A data frame reporting `event`, `mutant`, `fired`, and `time` for each
event type, for each species, and for each sampled time. Whenever, the
simulation has epigenetic states, the data frame also contains the
column `epistate`.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.1, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# place an "A[E1]" cell in the tissue
sim$place_cell("A[E1]", 500, 500)

# set delta time between species counting to 30
sim$history_delta <- 30

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# get the number of event fired per event and species
sim$get_firing_history()
#>           event mutant epistate fired     time
#> 1        deaths      A       E1    33 30.00000
#> 2  duplications      A       E1    62 30.00000
#> 3      switches      A       E1     3 30.00000
#> 4        deaths      A       E2     1 30.00000
#> 5  duplications      A       E2     2 30.00000
#> 6      switches      A       E2     0 30.00000
#> 7        deaths      A       E1   406 60.00000
#> 8  duplications      A       E1   674 60.00000
#> 9      switches      A       E1    31 60.00000
#> 10       deaths      A       E2    14 60.00000
#> 11 duplications      A       E2    90 60.00000
#> 12     switches      A       E2    13 60.00000
#> 13       deaths      A       E1   717 70.00099
#> 14 duplications      A       E1  1166 70.00099
#> 15     switches      A       E1    54 70.00099
#> 16       deaths      A       E2    21 70.00099
#> 17 duplications      A       E2   160 70.00099
#> 18     switches      A       E2    22 70.00099
```
