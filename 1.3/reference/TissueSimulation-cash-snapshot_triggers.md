# The snapshot triggers

This property is the tissue simulation snapshot trigger list

## Details

This property is a named list containing three values at most:
`time interval`, `clock interval`, and `number of cells`. They
represents the maximum computation time, the maximum simulation time,
and the maximum difference in the number of tumour cells between two
snapshots, respectively.

Notice that this property differs from
[`TissueSimulation$history_delta`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-history_delta.md).

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_code.md),
[`TissueSimulation$get_snapshot_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_snapshot_info.md),
[`recover_simulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/recover_simulation.md),
[`TissueSimulation$history_delta`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-history_delta.md)

## Examples

``` r
# set the random seed
set.seed(0)

# build a simulation
sim <- TissueSimulation()

# add a mutant
sim$add_mutant("A", list(duplication = 3, death = 1))

# place a cell
sim$place_cell("A", 500, 500)

# get snapshot information
sim$get_snapshot_info()
#>                  time clock cells
#> 1 2026-08-14 00:45:27     0     1
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d2e3a9.dat

# get the simulation's snapshot triggers
sim$snapshot_triggers
#> list()

# take a new snapshot every 1000 new tumour cells
sim$snapshot_triggers <- list("number of cells" = 1000)

# get new simulation's snapshot triggers
sim$snapshot_triggers
#> $`number of cells`
#> [1] 1000
#> 

# let the simulation evolve until consists of 5000 cells
sim$run_up_to_size("A", 4000)
#> 
 [███████████-----------------------------] 25% [00m:00s] Saving snapshot       

 [█████████████████████-------------------] 50% [00m:00s] Saving snapshot       

 [███████████████████████████████---------] 75% [00m:00s] Saving snapshot       

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# get snapshot information
sim$get_snapshot_info()
#>                  time    clock cells
#> 1 2026-08-14 00:45:27 0.000000     1
#> 2 2026-08-14 00:45:27 4.858114  1001
#> 3 2026-08-14 00:45:27 6.086508  2001
#> 4 2026-08-14 00:45:27 6.975337  3001
#> 5 2026-08-14 00:45:27 7.718603  4000
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d2e3a9.dat
#> 2 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d32334.dat
#> 3 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d37ed0.dat
#> 4 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d3f2ca.dat
#> 5 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d477d6.dat

# set new snapshot triggers: get a snapshot every 30 seconds, every 10
# simulated time units, and every 1000 new cells (because set previously)
sim$snapshot_triggers <- list("time interval" = as.difftime(30, units = "secs"),
                              "clock interval" = 10)

# get current simulation's snapshot triggers
sim$snapshot_triggers
#> $`time interval`
#> Time difference of 30 secs
#> 
#> $`clock interval`
#> [1] 10
#> 
#> $`number of cells`
#> [1] 1000
#> 

# take a snapshot every 5 simulated time units and reset the other
# trigger types
sim$snapshot_triggers <- list("time interval" = NULL,
                              "clock interval" = 5,
                              "number of cells" = NULL)

# get new simulation's snapshot triggers
sim$snapshot_triggers
#> $`clock interval`
#> [1] 5
#> 

# let the simulation evolve for other 20 time units
sim$run_up_to_time(sim$get_clock() + 20)
#> 
 [███████████████████---------------------] 45% [00m:00s] Saving snapshot       

 [█████████████████████-------------------] 52% [00m:00s] Cells: 19263          

 [█████████████████████████---------------] 62% [00m:01s] Cells: 29168          

 [██████████████████████████--------------] 63% [00m:01s] Saving snapshot       

 [████████████████████████████------------] 69% [00m:02s] Cells: 37494          

 [██████████████████████████████----------] 74% [00m:03s] Cells: 44851          

 [████████████████████████████████--------] 79% [00m:04s] Cells: 51024          

 [█████████████████████████████████-------] 81% [00m:05s] Saving snapshot       

 [██████████████████████████████████------] 83% [00m:05s] Cells: 57142          

 [███████████████████████████████████-----] 87% [00m:06s] Cells: 62613          

 [█████████████████████████████████████---] 90% [00m:07s] Cells: 67793          

 [█████████████████████████████████████---] 92% [00m:08s] Cells: 72397          

 [██████████████████████████████████████--] 94% [00m:09s] Cells: 75258          

 [███████████████████████████████████████-] 97% [00m:10s] Cells: 79675          

 [████████████████████████████████████████] 99% [00m:11s] Cells: 83296          

 [████████████████████████████████████████] 100% [00m:11s] Saving snapshot      


# get snapshot information
sim$get_snapshot_info()
#>                  time     clock cells
#> 1 2026-08-14 00:45:27  0.000000     1
#> 2 2026-08-14 00:45:27  4.858114  1001
#> 3 2026-08-14 00:45:27  6.086508  2001
#> 4 2026-08-14 00:45:27  6.975337  3001
#> 5 2026-08-14 00:45:27  7.718603  4000
#> 6 2026-08-14 00:45:27 12.718639 14173
#> 7 2026-08-14 00:45:29 17.718647 30895
#> 8 2026-08-14 00:45:32 22.718651 54792
#> 9 2026-08-14 00:45:39 27.718607 85008
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d2e3a9.dat
#> 2 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d32334.dat
#> 3 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d37ed0.dat
#> 4 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d3f2ca.dat
#> 5 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575d477d6.dat
#> 6 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575dd6ff8.dat
#> 7 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f575f642b9.dat
#> 8 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f57629e5a3.dat
#> 9 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_e780bb93/snapshot_658f57691bd73.dat
```
