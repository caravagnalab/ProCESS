# Getting the sampled cells' CNAs

This method returns the CNAs of the sample cells.

## Value

A data frame reporting `cell_id`, `type` (`"A"` for amplifications and
`"D"` for deletions), `chr`, `begin` (i.e., the first CNA locus in the
chromosome), `end` (i.e., last CNA locus in the chromosome), `allele`,
`src allele` (the allele origin for amplifications, `NA` for deletions),
and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
`"pre-neoplastic"`).

## Details

This method builds a data frame representing all the CNAs in the cells
sampled during the simulation and represented by the leaves of the
phylogenetic forest.

## See also

[`PhylogeneticForest$get_sampled_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md),
[`PhylogeneticForest$get_cell_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_cell_CNAs.md),
[`PhylogeneticForest$get_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_cell_mutations.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# build a tissue simulation with epigenetic states "E1" and "E2"
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01, E1 = 0.01)))

# place a cell of species "A[E1]" in position (500,500)
sim$place_cell("A[E1]", 500, 500)

# run the simulation until "A[E2]" accounts for less than 1000 cells
sim$run_up_to_size("A[E2]", 1000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# sample the tissue
sim$sample_cells("Sample_1", c(475, 475), c(525, 525))

# build the sample forest
sample_forest <- sim$get_sample_forest()

# initialize a mutation engine with the "demo" setup
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                                                        

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                                                                       

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                                                             

 [█████████████---------------------------] 30% [00m:01s] Loading RS index                                                                                                            

 [█████████████████████████---------------] 61% [00m:02s] Loading RS index                                                                                                            

 [████████████████████████████████████----] 89% [00m:03s] Loading RS index                                                                                                            

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                                                            

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                                                             

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                                                            


# add the genomic characterisation for the mutant "A"
m_engine$add_mutant("A",
                    list("E1" = c(SNV = 1e-7, indel = 1e-8),
                         "E2" = c(SNV = 3e-7, CNA = 1e-11)),
                    list(SNV("22", 10510210, "C", allele = 1),
                         CNA("D", "22", 5010000, 200000,
                             allele = 1)))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                                                                                          

 [████████████████████████████████████████] 100% [00m:00s] "A"'s SIDs validated                                                                                                       


# add the exposure
m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))

# build the phylogenetic forest
phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                                                                                            

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                                                                                           


# get the CNAs in the sampled cells
CNAs <- phylo_forest$get_sampled_cell_CNAs()

# show the first CNAs in the sampled cells
CNAs %>% head()
#>   cell_id type chr   begin     end allele src.allele  class
#> 1   19170    D  22 5010000 5209999      1         NA driver
#> 2   11472    D  22 5010000 5209999      1         NA driver
#> 3   17632    D  22 5010000 5209999      1         NA driver
#> 4   17633    D  22 5010000 5209999      1         NA driver
#> 5   20890    D  22 5010000 5209999      1         NA driver
#> 6   23064    D  22 5010000 5209999      1         NA driver
```
