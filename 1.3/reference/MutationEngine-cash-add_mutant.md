# Adding a mutant specification

This method adds a mutant specification to the mutation engine.

## Arguments

- mutant_name:

  The mutant name.

- passenger_rates:

  The list of the passenger rates whose names are the epigenetic states
  of the species or a single rate, if the mutant does not have an
  epigenetic state.

- drivers:

  The list of the driver SNVs, indels, CNAs, and the whole genome
  doubling events (WGD) characterizing the mutant (optional).

## Details

The users must use it to specify the name and the genomic
characterisation (i.e., SNVs, indels, CNAs, and whole genome doubling
events (WGD)) of all the simulated mutants together with the mutation
rates of its species. The driver mutations are applied to the mutant
progenitor's genome respecting the specification order.

## See also

[`MutationEngine$change_rates_from()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-change_rates_from.md),
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# create a demonstrative mutation engine
m_engine <- MutationEngine(setup_code = "demo", quiet = TRUE)

# define a list of mutations
d_mutations <- list("DGCR8 P26L",
                    Mutation("22", 16085675, "GCCTCCCGA", "G"),
                    "EP300 S2346del",
                    WGD,
                    CNA(type = "A", chr = "22", from = 10303470,
                        len = 200000),
                    SNV("22", 23657587, "C"),
                    CNA("D", "22", 5010000, 200000))

# add the mutant "A" characterized by the mutations in `d_mutations`. The
# mutations are applied according to `d_mutations`'s order. The mutant has
# one epigenetic states and its species "A[E1]" and "A[E2]" have passenger
# SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
# 1e-11, respectively.
m_engine$add_mutant("A", passenger_rates =c(SNV = 1e-9, indel = 1e-10),
                    drivers = d_mutations)
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                                                  

 [█---------------------------------------] 0% [00m:00s] Found 22                                                                             

 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                           

 [█████████████---------------------------] 32% [00m:01s] Reading 22                                                                          

 [███████████████████████████-------------] 65% [00m:02s] Reading 22                                                                          

 [███████████████████████████████████████-] 97% [00m:03s] Reading 22                                                                          

 [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated                                                               


m_engine
#> MutationEngine
#>  Passenger rates
#>    "A":
#>       [0,inf): {SNV: 1e-09, indel: 1e-10}
#> 
#>  Driver mutations
#>    "A":
#>        DGCR8 P26L (chr22(20073563)[C>T]) on random allele
#>        (chr22(16085675)[GCCTCCCGA>G]) on random allele
#>        EP300 S2346del (chr22(41574750)[TTC>]) on random allele
#>        Whole genome duplication
#>        CNA("A",chr22(10303470), len: 200000)
#>        (chr22(23657587)[G>C]) on random allele
#>        CNA("D",chr22(5010000), len: 200000)
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# add the mutant "B" characterized by the mutations in `d_mutations`. The
# mutations are applied according to `d_mutations`'s order. The mutant has
# two epigenetic states and its species "B[E1]" and "B[E2]" have passenger
# SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
# 1e-11, respectively.
m_engine$add_mutant("B", list("E1" = c(SNV = 1e-9, indel = 1e-10),
                              "E2" = c(SNV = 3e-8, CNA = 1e-11)),
                    drivers = d_mutations)
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs                                                                  

 [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated                                                               


m_engine
#> MutationEngine
#>  Passenger rates
#>    "A":
#>       [0,inf): {SNV: 1e-09, indel: 1e-10}
#>    "B[E1]":
#>       [0,inf): {SNV: 1e-09, indel: 1e-10}
#>    "B[E2]":
#>       [0,inf): {SNV: 3e-08, CNA: 1e-11}
#> 
#>  Driver mutations
#>    "A":
#>        DGCR8 P26L (chr22(20073563)[C>T]) on random allele
#>        (chr22(16085675)[GCCTCCCGA>G]) on random allele
#>        EP300 S2346del (chr22(41574750)[TTC>]) on random allele
#>        Whole genome duplication
#>        CNA("A",chr22(10303470), len: 200000)
#>        (chr22(23657587)[G>C]) on random allele
#>        CNA("D",chr22(5010000), len: 200000)
#>    "B":
#>        DGCR8 P26L (chr22(20073563)[C>T]) on random allele
#>        (chr22(16085675)[GCCTCCCGA>G]) on random allele
#>        EP300 S2346del (chr22(41574750)[TTC>]) on random allele
#>        Whole genome duplication
#>        CNA("A",chr22(10303470), len: 200000)
#>        (chr22(23657587)[G>C]) on random allele
#>        CNA("D",chr22(5010000), len: 200000)
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 
```
