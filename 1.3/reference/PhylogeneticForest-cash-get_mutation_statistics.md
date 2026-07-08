# Getting the statistics about mutations on each node

This method returns a dataframe reporting the statistics about mutations
on each node.

## Value

A dataframe consisting of five columns `cell_id`, `new_SIDs`,
`new_CNAs`, `total_SIDs`, and `total_CNAs`. Each row represents a node
in the phylogenetic forest and reports the identifier of the
corresponding cell and contains the number of mutations (`new_SIDs`) and
CNAs (`new_CNAs`) appearing for the first time on the cell. Moreover, it
show the total number of mutations and CNAs on the cell (`total_SIDs`
and `total_CNAs`, respectively).
