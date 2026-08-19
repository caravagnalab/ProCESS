# Building a frame generator representing the tissue configuration

This function builds a frame generator for
[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)
to represent each snapshot of `simulation` as its tissue configuration.

## Usage

``` r
get_tissue_frame_gen(simulation)
```

## Arguments

- simulation:

  The tissue simulation for which the frame generator is built.

## Value

A named list of two functions: `frame_generator` and `cleanup_function`.
The function `frame_generator` is a frame generator for
[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md)
that represents each snapshot of `simulation` as both its tissue
configuration. The function `cleanup_function` will cleanup the
temporary files required to execute `frame_generator` and must be called
after the frame generation has been completed.

## See also

[`build_snapshot_video()`](https://caravagnalab.github.io/ProCESS/1.3/reference/build_snapshot_video.md),
[`get_tissue_forest_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_tissue_forest_frame_gen.md),
[`get_forest_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_forest_frame_gen.md)
