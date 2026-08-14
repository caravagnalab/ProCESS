# Getting snapshot information

This method returns the data frame of the snapshots

## Value

A data frame consisting of four columns: `time`, `clock`, `cells`, and
`file`. Each row represents a snapshot. The column `time` stores the
snapshot time. The columns `clock` and `cells` maintain the simulated
time and the number of tumour cells at the snapshot time, respectively.
Finally, the column `file` contains the snapshot file path.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_code.md),
[`TissueSimulation$snapshot_triggers`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-snapshot_triggers.md),
[`recover_simulation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/recover_simulation.md)

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
#> 1 2026-08-14 00:44:56     0     1
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_27549a7c/snapshot_658f573fa0298.dat

# take a new snapshot every 10 simulated time units
sim$snapshot_triggers <- list("clock interval" = 10)

# let the simulation evolve up to time 30
sim$run_up_to_time(30)
#> 
 [██████████████--------------------------] 33% [00m:00s] Saving snapshot       

 [███████████████████---------------------] 46% [00m:00s] Cells: 18096          

 [███████████████████████-----------------] 56% [00m:01s] Cells: 28010          

 [██████████████████████████--------------] 63% [00m:02s] Cells: 36544          

 [███████████████████████████-------------] 66% [00m:03s] Saving snapshot       

 [████████████████████████████------------] 68% [00m:03s] Cells: 43803          

 [█████████████████████████████-----------] 72% [00m:04s] Cells: 50272          

 [███████████████████████████████---------] 76% [00m:05s] Cells: 56415          

 [█████████████████████████████████-------] 80% [00m:06s] Cells: 62001          

 [██████████████████████████████████------] 83% [00m:07s] Cells: 67129          

 [███████████████████████████████████-----] 85% [00m:08s] Cells: 72087          

 [███████████████████████████████████-----] 87% [00m:09s] Cells: 76144          

 [█████████████████████████████████████---] 90% [00m:10s] Cells: 80809          

 [█████████████████████████████████████---] 92% [00m:11s] Cells: 85139          

 [██████████████████████████████████████--] 94% [00m:12s] Cells: 89160          

 [███████████████████████████████████████-] 96% [00m:13s] Cells: 93085          

 [████████████████████████████████████████] 98% [00m:14s] Cells: 96893          

 [████████████████████████████████████████] 99% [00m:15s] Cells: 100741         

 [████████████████████████████████████████] 100% [00m:15s] Saving snapshot      


# get snapshot information
sim$get_snapshot_info()
#>                  time    clock  cells
#> 1 2026-08-14 00:44:56  0.00000      1
#> 2 2026-08-14 00:44:56 10.00003   7879
#> 3 2026-08-14 00:44:59 20.00003  41040
#> 4 2026-08-14 00:45:12 30.00000 101378
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_27549a7c/snapshot_658f573fa0298.dat
#> 2 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_27549a7c/snapshot_658f573fe3d26.dat
#> 3 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_27549a7c/snapshot_658f574311acc.dat
#> 4 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_27549a7c/snapshot_658f574f14ba2.dat
```
