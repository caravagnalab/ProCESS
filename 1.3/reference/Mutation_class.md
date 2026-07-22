# Either an SNV or an indel

A class to represents SNVs and indels

## Details

The objects of this class represent either an SNV or an indel and can be
created by
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)
and
[`Mutation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md).
They provide the following methods:

- [`get_alt()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_alt.md)
  returns the sequence as altered by the mutation.

- [`get_cause()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_cause.md)
  returns the cause of the mutation.

- [`get_chromosome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_chromosome.md)
  returns the chromosome in which the mutation occurred.

- [`get_dataframe()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_dataframe.md)
  returns a data frame representing the mutation.

- [`get_position_in_chromosome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_position_in_chromosome.md)
  returns the position in the chromosome of the mutation.

- [`get_ref()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation-cash-get_ref.md)
  returns the sequence as it was before the mutation.

## See also

[`Mutation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md)
