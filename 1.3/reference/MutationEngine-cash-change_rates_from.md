# Changing the passenger rates from a specified time

This method changes the passenger rates from a specified time.

## Arguments

- time:

  The time from which the passenger rates are set.

- mutant_name:

  The mutant name.

- passenger_rates:

  The list of the passenger rates whose names are the epigenetic states
  of the species or a single rate, if the mutant does not have an
  epigenetic state.

## Details

This method changes the passenger rates from a specified time. The rates
before the specified time and those of the unspecified epigenetic states
are not affected.

## See also

[`MutationEngine$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-add_mutant.md),
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# create a demonstrative mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 32% [00m:01s] Loading RS index                                                                    

 [██████████████████████████--------------] 64% [00m:02s] Loading RS index                                                                    

 [██████████████████████████████████████--] 94% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


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
# two epigenetic states and its species "A[E1]" and "A[E2]" have passenger
# SNV rates 1e-9 and 3e-8, respectively, and passenger CNA rates 0 and
# 1e-11, respectively.
m_engine$add_mutant("A", list("E1" = c(SNV = 1e-9, indel = 1e-10),
                              "E2" = c(SNV = 3e-8, CNA = 1e-11)),
                    drivers = d_mutations)
#> 
 [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs                                                                  

 [████████████████████████████████████████] 100% [00m:00s] "A"'s SIDs validated                                                               


m_engine
#> MutationEngine
#>  Passenger rates
#>    "A[E1]":
#>       [0,inf): {SNV: 1e-09, indel: 1e-10}
#>    "A[E2]":
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
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# change the rates of "A[E1]" from time 10
m_engine$change_rates_from(10, "A", list("E1" = c(SNV = 2e-9, indel = 4e-9)))

# ... and those of "A[E2]" from time 13
m_engine$change_rates_from(13, "A", list("E2" = c(SNV = 2e-9)))

m_engine
#> MutationEngine
#>  Passenger rates
#>    "A[E1]":
#>       [0,10): {SNV: 1e-09, indel: 1e-10}
#>       [10,inf): {SNV: 2e-09, indel: 4e-09}
#>    "A[E2]":
#>       [0,13): {SNV: 3e-08, CNA: 1e-11}
#>       [13,inf): {SNV: 2e-09}
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
```
