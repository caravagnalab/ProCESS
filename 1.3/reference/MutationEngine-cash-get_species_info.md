# Getting the registered species and their rates

This method returns the registered species and their rates.

## Value

A data frame containing the registered species rates.

## Details

The registered species and their rates during the simulation are
returned in a data frame. The column `mutant` contains the mutant names;
the columns `time`, `SNV_rate`, `indel_rate`, and `CNA_rate` store the
time from which rates hold, and the corresponding the SNV, indel, and
CNA rates, respectively. If the simulation has epigenetic states, then
the data frame also contains the column `epistate` to store the species
epigenetic state names.

## See also

[`SampleForest$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_species_info.md),
[`PhylogeneticForest$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_species_info.md),
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# get the active germline subject data frame
head(m_engine$get_species_info(), 5)
#> [1] mutant     time       SNV_rate   indel_rate CNA_rate  
#> <0 rows> (or 0-length row.names)
```
