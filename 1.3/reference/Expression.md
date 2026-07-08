# Represent an expression of value and variable

The objects of this class are polynomial expression of variables and
values built by using the binary operators `+`, `-`, and `*` with the
usual semantics.

An expression is one of the following object:

- a variable;

- a numeric value, e.g., `3.4`;

- the sum of two expressions;

- the subtraction of two expressions;

- the multiplication of two expressions.

## See also

[`Variable()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Variable.md)

## Examples

``` r
# build a simulation and add two species to it
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.2, death = 0.1))

# build an expression
sim$var("A") - 2 * sim$var("Time") * sim$var("A.duplications") + 3.4
#> |A|-2*Time*|A.duplications|+3.4

# R variables storing logic variables can also be used in expressions
v_time <- sim$var("Time")
sim$var("A") - 2 * v_time * sim$var("A.duplications") + 3.4
#> |A|-2*Time*|A.duplications|+3.4

# the logic expression can be stored in an R variable
v_exp <- sim$var("A") - 2 * v_time * sim$var("A.duplications") + 3.4
v_exp
#> |A|-2*Time*|A.duplications|+3.4
```
