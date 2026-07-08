# Set the rate of an event

This method sets the species' rate of an event.

## Arguments

- species:

  The species whose rate must be set.

- event_name:

  The name of the event whose rate should be set.

- dest:

  Either the species or the epigenetic state of one of the children due
  to the event (to be specified for the epigenetic switch event only).

- rate:

  The rate of the event.

## See also

[`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md),
[`TissueSimulation$get_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_rates.md)

## Examples

``` r
# create a simulation with epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2", "E3"))

# add mutant "A"
sim$add_mutant("A")

# set the duplication and death rates of the species "A[E1]"
sim$set_rate("A[E1]", "duplication", 0.1)
sim$set_rate("A[E1]", "death", 0.2)

# setting the switch event rates from "A[E1]" to "A[E2]" and "A[E3]"
sim$set_rate("A[E1]", "switch", "E2", 0.0001)
sim$set_rate("A[E1]", "switch", "A[E3]", 0.0002)

sim$get_rates()
#>   mutant epistate       event first.child.epistate  rate
#> 1      A       E1       death                 <NA> 2e-01
#> 2      A       E1 duplication                   E1 1e-01
#> 3      A       E1      switch                   E2 1e-04
#> 4      A       E1      switch                   E3 2e-04
```
