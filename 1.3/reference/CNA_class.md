# A copy number alteration

A class to represent CNA

## Details

This class represents copy number alterations (CNA). The objects of this
class are built by
[`CNA()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md),
[`Amplification()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Amplification.md),
and
[`Deletion()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Deletion.md).
They provide the following methods:

- [`get_allele()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA-cash-get_allele.md)
  returns the identifier of the CNA allele.

- [`get_dataframe()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA-cash-get_dataframe.md)
  returns a data frame representing the CNA.

- [`get_length()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA-cash-get_length.md)
  returns the length of the CNA.

- [`get_position_in_chromosome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA-cash-get_length.md)
  returns the position of the CNA in its chromosome.

- [`get_src_allele()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA-cash-get_src_allele.md)
  returns the identifier of the allele from which the new allele is
  copied when the CNA is an amplification.

## See also

[`CNA()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md),
[`Amplification()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Amplification.md),
[`Deletion()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Deletion.md)
