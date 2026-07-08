# Death activation level

The number of cells that activates cell death in a species.

## Details

This value is the minimum number of cells that enables cell death in a
species. The cell of a species \$S\$ can die if and only if that \$S\$
has reached the death activation level at least once during the
simulation.

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

# get the simulation death activation level
sim$death_activation_level
#> [1] 1

# set the death activation level to 50
sim$death_activation_level <- 50
```
