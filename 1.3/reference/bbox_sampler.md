# Bounding box sampler

This function searches for a tissue rectangle that contains a given
number of cells of a certain type. The search is performed by sampling
the tissue uniformly until the constraint on the number of cells
contained in the rectangle is satisfied or the sampling has been
repeated for a specified number of times. When the constraint cannot be
satisfied, the rectangle having the maximal number of cells of the
specified type among those sampled is returned.

## Usage

``` r
bbox_sampler(simulation, which, n, n_w, n_h, nattempts = 100)
```

## Arguments

- simulation:

  A simulation object.

- which:

  The species name.

- n:

  The desired number of cells from species `which`.

- n_w:

  Width of the box.

- n_h:

  Height of the box.

- nattempts:

  The maximum number of attempts.

## Value

coordinates for a bounding box.

## Examples

``` r
sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.08, death = 0.01))
sim$place_cell("A", 500, 500)
sim$run_up_to_size("A", 25000)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    

bbox <- bbox_sampler(sim, "A", n = 2500, n_w = 50, n_h = 50)
sim$sample_cells("A", bbox$p, bbox$q)
plot_tissue(sim)
```
