# Retrieving the rates update history

This method retrieves the simulation rates update history.

## Value

A data frame containing the set simulation rates. If the simulation has
epigenetic states, the data frame has 6 columns: `time`, `mutant`,
`epistate`, `event`, `first.child.epistate`, and `rate`. The column
`time` contains the rate setting time. The columns `mutant` and
`epistate` store the mutant and the epigenetic state of the cell from
which the event may occur. The columns `event` and `rate` maintain the
name and the rate of the event. Finally, the column
`first.child.epistate` reports the epigenetic state of potential first
child due to the event. For instance, when the event is `duplication`,
the first child has the same epigenetic state of its parent. Instead,
when the event is `switch`, the column `first.child.epistate` contains
an epigenetic state diffent from that of the origin cell. In case of the
event `death`, the column `first.child.epistate` is set to NA. When the
simulation has no epigenetic states, the returned data frame exclusively
contains the columns `time`, `mutant`, `event`, and `rate`.

## See also

[`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md),
[`TissueSimulation$get_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_rates.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.1))

# place a cell of "A"
sim$place_cell("A", 500, 500)

# set the death activation level
sim$death_activation_level <- 100

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# set the death rate of "A" to 0.9
sim$set_rate("A", "death", 0.9)

# we changed our mind *before* running the simulation
# and we reset the death rate to 0.05
sim$set_rate("A", "death", 0.05)

# simulate up to time 80
sim$run_up_to_time(80)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# set the death rate to 0.5
sim$set_rate("A", "death", 0.5)

# simulate up to time 80+1
sim$run_up_to_time(sim$get_clock()+1)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# get the rates update history
sim$get_rates_update_history()
#>       time mutant       event rate
#> 1  0.00000      A       death 0.10
#> 2  0.00000      A duplication 0.20
#> 3 70.00092      A       death 0.05
#> 4 80.00073      A       death 0.50
```
