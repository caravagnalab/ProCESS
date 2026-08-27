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
#> 1 2026-08-26 22:27:18     0     1
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_6b30a416/snapshot_659f90bbb4c88.dat

# take a new snapshot every 10 simulated time units
sim$snapshot_triggers <- list("clock interval" = 10)

# let the simulation evolve up to time 30
sim$run_up_to_time(30)
#> 
 [██████████████--------------------------] 33% [00m:00s] Saving snapshot                                                                     

 [███████████████████---------------------] 46% [00m:00s] Cells: 18120                                                                        

 [███████████████████████-----------------] 56% [00m:01s] Cells: 28339                                                                        

 [██████████████████████████--------------] 63% [00m:02s] Cells: 36740                                                                        

 [███████████████████████████-------------] 66% [00m:03s] Saving snapshot                                                                     

 [████████████████████████████------------] 68% [00m:03s] Cells: 44060                                                                        

 [██████████████████████████████----------] 73% [00m:04s] Cells: 50739                                                                        

 [███████████████████████████████---------] 77% [00m:05s] Cells: 56964                                                                        

 [█████████████████████████████████-------] 80% [00m:06s] Cells: 62650                                                                        

 [██████████████████████████████████------] 83% [00m:07s] Cells: 68041                                                                        

 [███████████████████████████████████-----] 86% [00m:08s] Cells: 72957                                                                        

 [████████████████████████████████████----] 88% [00m:09s] Cells: 76572                                                                        

 [████████████████████████████████████----] 89% [00m:10s] Cells: 79554                                                                        

 [█████████████████████████████████████---] 91% [00m:11s] Cells: 83584                                                                        

 [██████████████████████████████████████--] 93% [00m:12s] Cells: 87614                                                                        

 [███████████████████████████████████████-] 95% [00m:13s] Cells: 91631                                                                        

 [███████████████████████████████████████-] 97% [00m:14s] Cells: 95376                                                                        

 [████████████████████████████████████████] 99% [00m:15s] Cells: 99083                                                                        

 [████████████████████████████████████████] 100% [00m:16s] Saving snapshot                                                                    


# get snapshot information
sim$get_snapshot_info()
#>                  time    clock  cells
#> 1 2026-08-26 22:27:18  0.00000      1
#> 2 2026-08-26 22:27:19 10.00003   7879
#> 3 2026-08-26 22:27:22 20.00003  41040
#> 4 2026-08-26 22:27:35 30.00000 101378
#>                                                                                                   file
#> 1 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_6b30a416/snapshot_659f90bbb4c88.dat
#> 2 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_6b30a416/snapshot_659f90bbf8a6f.dat
#> 3 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_6b30a416/snapshot_659f90bf1c2f8.dat
#> 4 /var/folders/tb/jqmdpgxs2t5129bny6pb96680000gn/T/ProCESS_alberto_6b30a416/snapshot_659f90cb8edfc.dat
```
