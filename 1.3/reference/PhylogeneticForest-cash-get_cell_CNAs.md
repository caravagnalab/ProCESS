# Getting the CNAs of a cell in the samples' phylogeny

This method returns the CNAs of one cell represented in the forest

## Arguments

- cell_id:

  The identifier of the cell whose CNAs are aimed.

## Value

A data frame reporting `cell_id`, `type` (`"A"` for amplifications and
`"D"` for deletions), `chr`, `begin` (i.e., the first CNA locus in the
chromosome), `end` (i.e., last CNA locus in the chromosome), `allele`,
`src allele` (the allele origin for amplifications, `NA` for deletions),
and `class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
`"pre-neoplastic"`).

## Details

This method builds a data frame representing all the CNAs in the cell
represented by nodes of the phylogenetic forest being either one of the
sampled cells or one of their ancestors.

## See also

[`PhylogeneticForest$get_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_cell_mutations.md),
[`PhylogeneticForest$get_sampled_cell_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_CNAs.md),
[`PhylogeneticForest$get_sampled_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md)

## Examples

``` r
# set the random seed for repeatability
set.seed(0)

# build a tissue simulation with epigenetic states "E1" and "E2"
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add mutant "A" and set its species rates
sim$add_mutant("A",
               list(E1 = list(duplication = 0.2, death = 0.05, E2 = 0.01),
                    E2 = list(duplication = 0.08, death = 0.01)))

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
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
#> Downloading signature files...
#> Signature file downloaded
#> Downloading driver mutation file...
#> Driver mutation file downloaded
#> Decompressing driver mutation file...done
#> Downloading passenger CNAs file...
#> Passenger CNAs file downloaded
#> Decompressing passenger CNAs file...done
#> Downloading germline...
#> Germline downloaded
#> Decompressing mutations...
#> done
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                   

 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22                                  

 [██████████████████████████████----------] 73% [00m:01s] Processing chr. 22                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                

#> done
#> Building repeated sequence index...
#> 
 [█---------------------------------------] 0% [00m:00s] Reading 22                                           

 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] Reading 22                                         

#> 
 [████████████████████████████████████████] 100% [00m:01s] RS index built                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Saving RS index                                      

 [█---------------------------------------] 0% [00m:00s] Saving RS index                                      

 [████████████████████████----------------] 58% [00m:02s] Saving RS index                                     
done
#> 
 [████████████████████████████████████████] 100% [00m:02s] RS index saved                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Saving germline                                      

 [████████████████████████████████████████] 100% [00m:00s] Germline saved                                     


# add the genomic characterisation for the mutant "A"
m_engine$add_mutant("A",
                    list("E1" = c(SNV = 1e-7, indel = 1e-8),
                         "E2" = c(SNV = 3e-7, CNA = 2e-8)),
                    list(SNV("22", 10510210, "C", allele = 1),
                         CNA("D", "22", 5010000, 200000,
                             allele = 1)))
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                  

 [█---------------------------------------] 0% [00m:00s] Found 22                                             

 [█---------------------------------------] 0% [00m:00s] Reading 22                                           

 [█████████████---------------------------] 32% [00m:01s] Reading 22                                          

 [██████████████████████████--------------] 63% [00m:02s] Reading 22                                          

 [███████████████████████████████████████-] 96% [00m:03s] Reading 22                                          

 [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated                               


# add the exposure
m_engine$add_exposure(c(ID1 = 1, SBS1 = 0.5, SBS2 = 0.5))

# build the phylogenetic forest
phylo_forest <- m_engine$place_mutations(sample_forest, 1, 1)
#> 
 [█---------------------------------------] 0% [00m:00s] Placing mutations                                    

 [████████████████████████████████████████] 100% [00m:00s] Mutations placed                                   


# load dplyr to use %>%
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# get the node corresponding to the youngest non-sampled cell
# belonging to "A[E2]"
node <- phylo_forest$get_nodes() %>%
  dplyr::filter(is.na(sample) & epistate == "E2") %>%
  dplyr::slice_max(birth_time, n = 1)

# gets the CNAs in it
phylo_forest$get_cell_CNAs(node$cell_id)
#>   cell_id type chr    begin      end allele src.allele     class
#> 1   22009    D  22  5010000  5209999      1         NA    driver
#> 2   22009    D  22 17550013 29065240      0         NA passenger
#> 3   22009    A  22 17900000 22735025      2          1 passenger
#> 4   22009    A  22 22387145 48423213      3          1 passenger
```
