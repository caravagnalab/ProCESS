# Getting the driver mutations

This method returns the applied driver mutations.

## Value

A data frame consisting in eight columns `mutant`, `order`, `type`,
`CNA_type`, `chr`, `start`, `end`, `ref`, `alt`, `allele`, `allele_srd`
and `code`. Each row in the data frame reports one driver mutations. The
fields `mutant` and `order` report the name of the mutant and the
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
