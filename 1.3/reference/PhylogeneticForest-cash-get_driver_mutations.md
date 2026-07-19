# Getting the driver mutations

This method returns the applied driver mutations.

## Value

A data frame consisting in eight columns `order`, `type`, `CNA_type`,
`chr`, `start`, `end`, `ref`, `alt`, `allele`, `allele_srd`, `cause`,
and `code`. Each row in the data frame reports one driver mutations. The
fields `cause` and `order` report the name of the mutant and the
application order among the mutant driver mutations, respectively. The
column `type` declares the mutation type and contains `SID`, `CNA`, or
`WGD` when the mutation is an SNV/indel, a CNA, or a whole genome
duplication, respectively. When the mutation is a CNA, the `CNA_type`
can either be `A` (i.e., amplification) or `D` (i.e., deletion). When
the mutation is not a WGD, the fields `chr`, `start`, and `end` contains
the mutation chromosome, the initial and the final position on the
chromosome, respectively. When the mutation is a SID , the fields `ref`
and `alt` contains the mutation reference genome and alternate
sequences, respectively. When the mutation is a SID or a CNA deletion,
the field `allele` stores the allele in which the mutation was applied.
When instead the mutation is a CNA amplification, the fields `allele`
and `src_allele` reports the identifiers of the new allele and of the
original allele, respectively. In all the remaining cases, the fields
contains the value `NA`. Finally, the column `code` reports the mutation
code.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# show information about samples
forest$get_driver_mutations()
#>   order type CNA_type  chr    start      end       ref  alt allele src.allele
#> 1     1  SID     <NA>   22 46510210 46510210         A    C      1         NA
#> 2     2  SID     <NA>   22 16085675 16085683 GCCTCCCGA    G      0         NA
#> 3     3  CNA        A   22 10303470 10503469      <NA> <NA>      2          0
#> 4     4  CNA        D   22  5010000  5209999      <NA> <NA>      1         NA
#> 5     5  WGD     <NA> <NA>       NA       NA      <NA> <NA>     NA         NA
#> 6     1  SID     <NA>   22 20073563 20073563         C    T      1         NA
#> 7     2  SID     <NA>   22 12028576 12028576         N    G      0         NA
#>   cause       code
#> 1     A       <NA>
#> 2     A       <NA>
#> 3     A       <NA>
#> 4     A       <NA>
#> 5     A       <NA>
#> 6     B DGCR8 P26L
#> 7     B       <NA>
```
