# Formula-based simulation constraints

> *Note:* This article presents advanced topics on tissue simulation.
> Refer to [this
> article](https://caravagnalab.github.io/ProCESS/1.3/articles/tissue_simulation.md)
> for an introduction on the subject.

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS implements a first-order unquantified logic having variables
representing the cardinality of the species, the number of events fired
in a species (being duplications, deaths, or switches), and the
simulation time. These variables and real values are summed by `+`,
subtracted by `-`, and multiplied by `*` and form expressions. The
expressions are then compared with the standard semantics by `>`, `>=`,
`==`, `!=`, `<=`, and `<` to form relations. A formula in this language
is either a relation, the conjunction of two formulas (`&`), or the
disjunction of two formulas (`|`).

Any formula in the above language expresses a condition on the
simulation status. It can be used as a parameter of the method
[`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)
to let the simulation evolve until the condition no longer holds.

### Variables

The variables represent one of the following quantities:

- the cardinality of a species;
- the number of fired deaths, duplications, or switches in a species;
- the elapse simulation time.

All the above variables can be built by using the method
[`TissueSimulation$var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md).
When the parameter is the string `"Time"`, the elapsed simulation time
variable is returned.

``` r

library("ProCESS")

# set the seed of the random number generator
set.seed(0)

# build a tissue simulation with two epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))

# add the mutant "A"
sim$add_mutant(name = "A", list(E1 = list(duplication = 0.2, death = 0.1,
                                          E2 = 0.01),
                                E2 = list(duplication = 0.08, death = 0.01,
                                          E1 = 0.01)))

# get the variable representing the simulation time
v_time <- sim$var("Time")

v_time
#> Time
```

When the parameter is the name of a species, a variable representing the
cardinality of the species is built.

``` r

# get the variable representing the cardinality of A[E1] in sim
va1 <- sim$var("A[E1]")
va1
#> |A[E1]|

# get the variable representing the cardinality of A[E2] in sim
va2 <- sim$var("A[E2]")
va2
#> |A[E2]|
```

Finally, when the parameter is the name of a species followed by a `.`
and the name of an event among `deaths`, `duplications`, or `switches`,
[`TissueSimulation$var()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-var.md)
returns the variable associated with the number of the corresponding
event in the specified species.

``` r

# get the variable representing the number of epigenetic
# switches from A[E1]
va1_s <- sim$var("A[E1].switches")
va1_s
#> |A[E1].switches|

# get the variable representing the number of duplications in A+
sim$var("A[E1].duplications")
#> |A[E1].duplications|

# get the variable representing the number of deaths in A+
sim$var("A[E1].deaths")
#> |A[E1].deaths|
```

### Expressions and Formulas

An expression is one of the following objects:

- a variable, e.g., `sim$var("A[E1]")`;
- a numeric value, e.g., `3.4`;
- the sum of two expressions, e.g., `sim$var("A[E1]") + 3.4`;
- the subtraction of two expressions, e.g., `sim$var("A[E1]") - 3.4`;
- the multiplication of two expressions, e.g., `sim$var("A[E1]") * 3.4`.

Two expression can be related by `<=`, `<`, `==`, `!=`, `>` and `>=`.

A formula is:

- a relation among two expressions, e.g., `sim$var("A[E1]")>=2`;
- the conjunction of two formulas, e.g.,
  `sim$var("A[E1]")>=2 & sim$var("A[E1]")<=500`;
- the disjunction of two formulas, e.g.,
  `sim$var("A[E1]")>=2 | sim$var("A[E1]")<=500`.

### The method [`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)

The method
[`TissueSimulation$run_until()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-run_until.md)
takes a formula as a parameter and runs the simulation until the formula
no longer holds.

We can build a condition stating that the cardinality of `A[E1]` doubles
that of `A[E2]`

``` r

c1 <- va1 >= 2 * va2
c1
#> |A[E1]|>=2*|A[E2]|
```

…, a condition that holds when there are more than 100000 live cells of
mutant `A`

``` r

c2 <- va1 + va2 > 1e5
c2
#> |A[E1]|+|A[E2]|>100000
```

…, a condition that constrains the simulated time

``` r

# build a condition that holds when 40 time unit have been
# simulated at least
c4 <- v_time >= 40
c4
#> Time>=40
```

…, or deal with the number of species events. For instance, we can build
a condition that holds when less than 4000 epigenetic switches from the
species `A[E1]` have occurred.

``` r

# build a condition that holds when less than 4000
# epigenetic switches from the species A[E1] have occurred
c3 <- va1_s < 4e3
c3
#> 4000>|A[E1].switches|
```

Moreover, we can logically combine different conditions.

``` r

# build a condition that holds when c4 and at least one
# among c1, c2, and c3 hold
c5 <- c4 & (c1 | c2 | c3)
c5
#> Time>=40 and (|A[E1]|>=2*|A[E2]| or |A[E1]|+|A[E2]|>100000 or 4000>|A[E1].switches|)
```

Any simulation can be run until a condition holds.

``` r

# place the initial cell
sim$place_cell("A[E1]", 500, 500)

# run the simulation while c5 does not hold
sim$run_until(c5)
#>  [████████████████████████████████████████] 100% [00m:00s] Saving snapshot

sim
#> ──  ProCESS   D   S   M  ProCESS_20260708-173741 ─────── ▣  [1000x1000]  ⏱ 40 ──
#> 
#> ── Species: 2, with epigenetics
#>    
#>    =======  ====  ====  ======  ========
#>    species   λ     δ    counts     %    
#>    =======  ====  ====  ======  ========
#>      A[E1]  0.20  0.10    84    74.33628
#>      A[E2]  0.08  0.01    29    25.66372
#>    =======  ====  ====  ======  ========
#> 
#> ── Epigenetic switches
#>    
#>    =======  ====  =====
#>    species   ε     dest
#>    =======  ====  =====
#>      A[E1]  0.01  A[E2]
#>      A[E2]  0.01  A[E1]
#>    =======  ====  =====
#> 
#> ── Firings: 288 total
#> 
#>  Species A[E1]:  78 (deaths), 171 (duplications) and  13 (switches)
#>  Species A[E2]:   2 (deaths),  21 (duplications) and   3 (switches)
#> ✖ The simulation has no samples yet!
sim$get_clock()
#> [1] 40.00892
```
