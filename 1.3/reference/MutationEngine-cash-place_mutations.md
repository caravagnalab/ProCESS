# Placing the mutations

This methods places mutations on a sample forest.

## Arguments

- sample_forest:

  A sample forest.

- num_of_pre_neoplastic_SNVs:

  The number of pre-neoplastic SNVs.

- pre_neoplastic_SNV_signature_name:

  The name of the SNV signature for the pre-neoplastic SNV generation
  (optional).

- num_of_pre_neoplastic_indels:

  The number of pre-neoplastic indels.

- pre_neoplastic_indel_signature_name:

  The name of the indel signature for the pre_neoplastic indel
  generation.

- seed:

  The seed for random number generator (optional).

## Value

A phylogenetic forest whose structure corresponds to `sample_forest`.

## Details

Each node of a sample forest is labelled by the mutations occurring in
the cell represented by the node itself and produces a phylogenetic
forest.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# build a tissue simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell of species "A" in position (500,500)
sim$place_cell("A", 500, 500)

# run the simulation until "A" accounts for less than 50000 cells
sim$run_up_to_size("A", 50000)
#> 
 [███████████████████████████████████████-] 96% [00m:00s] Cells: 48230          

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot      


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 475),
                 upper_corner = c(500, 550))

# build the sample forest
sample_forest <- sim$get_sample_forest()

# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# add the mutant "A" to the engine
m_engine$add_mutant("A", c(SNV = 3e-9), list(SNV("22", 12028576, "G")))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs    

 [████████████████████████████████████████] 100% [00m:00s] "A"'s SIDs validated 


# add the default set of SNV signature coefficients
m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7, ID2 = 0.3, ID21 = 0.5,
                        ID3 = 0.2))

# place the mutations on the sample forest assuming 1000 pre-neoplastic
# SNVs and 500 indels
phylogenetic_forest <- m_engine$place_mutations(sample_forest, 1000, 500)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations      

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed     


phylogenetic_forest
#> PhylogeneticForest
#>   # of trees: 1
#>   # of nodes: 13794
#>   # of leaves: 3854
#>   samples: {"S1"}
#> 
#>   # of emerged SNVs and indels: 5753
#>   # of emerged CNAs: 0
#> 
```
