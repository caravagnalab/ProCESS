# Set the tissue simulation rates

This method sets tissue simulation rates.

## Arguments

- rates:

  Either a list or a data frame of the rates to be set.

## Details

This method can set the rates of multiple species. It accepts as
parameters both a list of a data frame specifying the rates to be set.

- The list must be a named list whose names correspond to the species
  whose rates must be set. Each element in the list represents the rates
  to be set for the corresponding species and is itself a named list
  whose names must belong to the set containing `"death"`,
  `"duplication"`, and, if the simulation includes epigenetic states,
  the names of the known epigenetic states. The values of the element
  named `"duplication”` and `“death”` represent the new species’
  duplication and death rates, respectively. Instead, the values whose
  names are among the known epigenetic states indicate the new switch
  rate to the specified epigenetic state (see below for examples).

- The data frame must contain at least three columns: `mutant`, `event`,
  and `rate`. If the simulation includes epigenetic states, the data
  frame must also include the columns `epistate` and
  `first.child.epistate`. Each row in the data frame declares the new
  rate value for an event in a species. The columns `mutant` and `event`
  represent the species mutants and event types. The values in the
  former column must be known mutant names, while those in the latter
  must be among `“death”`, `“duplication”`, and, when the simulation
  includes epigenetic states, `“switch”`. The columns `epistate` and
  `first.child.epistate` denote the epigenetic state of the species
  whose rate is to be set and the epigenetic state of the first child
  due to the event. When the event is `“duplication”`, the
  `first.child.epistate` value should be equal to the value in the
  column `epistate`; when the event is `“switch”`, it must be a known
  epigenetic state different from that reported in the column epistate;
  when the event is `“death”`, it can be `NA`. Finally, the column
  `rate` contains the new rate values and must be numerical.

Notice that no new epigenetic state nor mutants will be added to the
simulation by this method. Any mention to a non-existant mutant or
epigenetic state ends the execution with an error.

## See also

[`TissueSimulation$set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md),
[`TissueSimulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# add mutants "A", "B", and "C"
sim$add_mutants(c("A", "B", "C"))

# no rate is reported because none has been set yet
sim$get_rates()
#> [1] mutant event  rate  
#> <0 rows> (or 0-length row.names)

# setting "A"'s duplication and death rates and "B"'s death rate
# All remaining rates are set to zero by default
sim$set_rates(list(A = list(duplication = 0.3, death = 0.2),
                   B = list(death = 0.2)))

# get_rates() only reports the rates that have been explicitly set
sim$get_rates()
#>   mutant       event rate
#> 1      A       death  0.2
#> 2      A duplication  0.3
#> 3      B       death  0.2

# adding epigenetic states
sim$add_epigenetic_states(c("E1", "E2", "E3"))

# species before adding epigenetic states do not exists anymore
# and no rates have been set for the new species
sim$get_rates()
#> [1] mutant               epistate             event               
#> [4] first.child.epistate rate                
#> <0 rows> (or 0-length row.names)

# set some rates for the new species
sim$set_rates(list("A[E1]" = list(duplication = 0.3, death = 0.2,
                                  E2 = 0.01, "A[E3]" = 0.04),
                   "B[E2]" = list(death = 0.2, E3 = 0.01)))

sim$get_rates()
#>   mutant epistate       event first.child.epistate rate
#> 1      A       E1       death                 <NA> 0.20
#> 2      A       E1 duplication                   E1 0.30
#> 3      A       E1      switch                   E2 0.01
#> 4      A       E1      switch                   E3 0.04
#> 5      B       E2       death                 <NA> 0.20
#> 6      B       E2      switch                   E3 0.01
# create a simulation
sim <- TissueSimulation()

# add mutants "A", "B", and "C"
sim$add_mutants(c("A", "B", "C"))

# build a data frame for setting "A"'s duplication and death
# rates and "B"'s death rate as well
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

# setting the rates by using the data frame
sim$set_rates(df_rates)

# get_rates() only reports the set rates.
sim$get_rates()
#>   mutant       event rate
#> 1      A       death  0.2
#> 2      A duplication  0.3
#> 3      B duplication  0.2

# adding epigenetic states
sim$add_epigenetic_states(c("E1", "E2", "E3"))

# now we need to add the "epistate" and `first.child.epistate`
# columns to the data frame
df_rates[["epistate"]] <- c("E1", "E1", "E2")
df_rates[["first.child.epistate"]] <- c("E1", NA, "E2")

# load dplyr to simplify the next part of the example
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# we also may set some switch rates
df_rates <- df_rates %>%
  add_row(mutant = "A", epistate = "E1", event = "switch",
          first.child.epistate = "E2", rate = 0.01) %>%
  add_row(mutant = "A", epistate = "E1", event = "switch",
          first.child.epistate = "E3", rate = 0.04) %>%
  add_row(mutant = "B", epistate = "E2", event = "switch",
          first.child.epistate = "E3", rate = 0.01)

df_rates
#>   mutant       event rate epistate first.child.epistate
#> 1      A duplication 0.30       E1                   E1
#> 2      A       death 0.20       E1                 <NA>
#> 3      B duplication 0.20       E2                   E2
#> 4      A      switch 0.01       E1                   E2
#> 5      A      switch 0.04       E1                   E3
#> 6      B      switch 0.01       E2                   E3

sim$set_rates(df_rates)

# get_rates() only reports the set rates. All remaining rates are
# set to zero by default
sim$get_rates()
#>   mutant epistate       event first.child.epistate rate
#> 1      A       E1       death                 <NA> 0.20
#> 2      A       E1 duplication                   E1 0.30
#> 3      A       E1      switch                   E2 0.01
#> 4      A       E1      switch                   E3 0.04
#> 5      B       E2 duplication                   E2 0.20
#> 6      B       E2      switch                   E3 0.01
```
