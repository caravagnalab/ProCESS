# Getting the simulation rates

This method returns the rates of the simulation.

## Arguments

- complete:

  A Boolean flag to get also the rates that have not been set (default:
  FALSE).

## Value

A data frame containing the simulation rates. If the simulation has
epigenetic states, the data frame has 5 columns: `mutant`, `epistate`,
`event`, `first.child.epistate`, and `rate`. The columns `mutant` and
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
contains the columns `mutant`, `event`, and `rate`.

## Details

This method returns a data frame containing the simulation rates. A rate
is not included in the returned data frame if and only if it was not set
during the system specification. In these cases, the rate is assumed to
be 0 by default.

## See also

[`TissueSimulation$set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# add a mutant
sim$add_mutant("A")

# set its duplication rate
sim$set_rate("A", "duplication", 0.1)

# get the rates that have been set
sim$get_rates()
#>   mutant       event rate
#> 1      A duplication  0.1

# get all simulation rates
sim$get_rates(TRUE)
#>   mutant       event rate
#> 1      A duplication  0.1
#> 2      A       death  0.0

# add epigenetic states
sim$add_epigenetic_states(c("E1", "E2"))

# set some of the rates of "A[E1]" and "A[E2]"
sim$set_rate("A[E1]", "duplication", 0.1)
sim$set_rate("A[E2]", "duplication", 0.1)
sim$set_rate("A[E1]", "death", 0.2)

# get the rates that have been set
sim$get_rates()
#>   mutant epistate       event first.child.epistate rate
#> 1      A       E1       death                 <NA>  0.2
#> 2      A       E1 duplication                   E1  0.1
#> 3      A       E2 duplication                   E2  0.1

# get all simulation rates
sim$get_rates(TRUE)
#>   mutant epistate       event first.child.epistate rate
#> 1      A       E1 duplication                   E1  0.1
#> 2      A       E1       death                 <NA>  0.2
#> 3      A       E1      switch                   E2  0.0
#> 4      A       E2 duplication                   E2  0.1
#> 5      A       E2       death                 <NA>  0.0
#> 6      A       E2      switch                   E1  0.0
```
