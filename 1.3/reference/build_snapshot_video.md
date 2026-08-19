# Building a video of the snapshots

This function builds a video of the snapshots.

## Usage

``` r
build_snapshot_video(simulation)
```

## Arguments

- simulation:

  A simulation object.

- output_file:

  The path of the output video. When it is set to `NULL`, the output
  video path has the format `<simulation path>_evolution.mp4` (default:
  `NULL`).

- width:

  The width of the video (default: `800`).

- height:

  The height of the video (default: `600`).

- framerate:

  The video framerate in frame/sec (default: `1`).

- res:

  The video resolution (default: `150`).

- frame_generator:

  The function used to plot each frame of the video. When it is set to
  `NULL`, this function uses the function `frame_generator` built by
  [`get_tissue_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_tissue_frame_gen.md)
  (default: `NULL`).

- pauses_on_event:

  A named list specifying the pauses on event in the output video. The
  names of the list must be among "mutant emerged" and "sampling". The
  values can either be a numeric value, a `difftime` object, or a named
  list. The numeric value represents the pause length in number of
  frames. The `difftime` object denotes the pause length. Finally, the
  named list described the pauses for specific events. When the name of
  the element is "mutant emerged", the names of the sub-list are among
  the simulated mutants. When, instead, the name of the element is
  "sampling", the names of the sub-list are among the collected samples.
  In both the cases, the values represent the pause lengths for the
  specific event and can be either a numeric value or a `difftime`
  object as for the generic specification. When `pauses_on_event` is
  `NULL`, a 3-seconds pause is added after any new mutant appearance and
  any sample collection (default: `NULL`).

- pauses_on_frame:

  A list specifying the pauses on frame in the output video. Each
  element of the list is a named list whose names are `frame` and
  `length`. The element `frame` is the video frame from which the pause
  is aimed. Instead, `length` is the pause's length expressed in either
  number of frames, when `length` is a numeric value, or a `difftime`.
  When `pauses_on_frame` is NULL, no pauses are added to the output
  video (default: `NULL`).

- quiet:

  A Boolean flag to enable/disable the messages (default: `FALSE`).

- workers:

  The number of parallel processes generating frames. This parameter is
  used only when the packages
  [furrr::furrr-package](https://furrr.futureverse.org/reference/furrr-package.html)
  and
  [progressr::progressr-package](https://progressr.futureverse.org/reference/progressr.html)
  are installed. When it is set to `NULL`, the function uses as many
  processes as the number of available processors minus one (default:
  `NULL`).

## Value

The name of the produced video file path.

## Details

This function builds a video of the snapshots. It is available only if
the package av::av-package is installed.

## See also

[`vignette("videos")`](https://caravagnalab.github.io/ProCESS/1.3/articles/videos.md),
[`plot_tissue()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md),
[`get_tissue_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_tissue_frame_gen.md),
[`get_forest_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_forest_frame_gen.md),
[`get_tissue_forest_frame_gen()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_tissue_forest_frame_gen.md)
