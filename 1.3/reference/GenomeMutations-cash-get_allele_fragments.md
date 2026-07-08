# Getting genome allele fragments

This method returns a data frame representing the allele fragments in
the genome.

## Value

A data frame consisting of 5 columns: `chr`, `allele`, `src_allele`,
`from`, and `size`. Each row represent a allele fragment. The columns
`chr`, and `allele` represent the fragment's chromosome and allele,
respectively. The column `allele_src` stores the allele from which the
allele of the fragment is derived. The columns `from` and `size`
maintain the first position of the fragment in the wild-type allele and
its size.
