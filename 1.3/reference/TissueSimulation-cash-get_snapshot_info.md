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
#> 1 2026-08-11 22:44:39     0     1
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_1b958e54/snapshot_658cb8a24df10.dat

# take a new snapshot every 10 simulated time units
sim$snapshot_triggers <- list("clock interval" = 10)

# let the simulation evolve up to time 30
sim$run_up_to_time(30)
#> 
 [██████████████--------------------------] 33% [00m:00s] Saving snapshot       

 [███████████████████---------------------] 46% [00m:00s] Cells: 18091          

 [███████████████████████-----------------] 56% [00m:01s] Cells: 28194          

 [██████████████████████████--------------] 63% [00m:02s] Cells: 36622          

 [███████████████████████████-------------] 66% [00m:03s] Saving snapshot       

 [████████████████████████████------------] 68% [00m:03s] Cells: 43951          

 [██████████████████████████████----------] 73% [00m:04s] Cells: 50536          

 [███████████████████████████████---------] 75% [00m:05s] Cells: 54248          

 [████████████████████████████████--------] 78% [00m:06s] Cells: 59760          

 [█████████████████████████████████-------] 81% [00m:07s] Cells: 65017          

 [██████████████████████████████████------] 84% [00m:08s] Cells: 69788          

 [███████████████████████████████████-----] 86% [00m:09s] Cells: 74295          

 [████████████████████████████████████----] 89% [00m:10s] Cells: 78727          

 [█████████████████████████████████████---] 91% [00m:11s] Cells: 82868          

 [██████████████████████████████████████--] 93% [00m:12s] Cells: 87125          

 [███████████████████████████████████████-] 95% [00m:13s] Cells: 90895          

 [███████████████████████████████████████-] 96% [00m:14s] Cells: 94527          

 [████████████████████████████████████████] 98% [00m:15s] Cells: 98337          

 [████████████████████████████████████████] 100% [00m:16s] Saving snapshot      


# get snapshot information
sim$get_snapshot_info()
#>                  time    clock  cells
#> 1 2026-08-11 22:44:39  0.00000      1
#> 2 2026-08-11 22:44:39 10.00003   7879
#> 3 2026-08-11 22:44:43 20.00003  41040
#> 4 2026-08-11 22:44:56 30.00000 101378
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_1b958e54/snapshot_658cb8a24df10.dat
#> 2 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_1b958e54/snapshot_658cb8a291c7d.dat
#> 3 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_1b958e54/snapshot_658cb8a5b9479.dat
#> 4 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_1b958e54/snapshot_658cb8b25a378.dat
```
