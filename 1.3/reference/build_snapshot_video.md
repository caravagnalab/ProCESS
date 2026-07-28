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

- plot_function:

  The function used to plot each frame of the video. When is set to
  `NULL`, the function
  [`plot_tissue()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)
  is used (default: `NULL`).

- width:

  The width of the video (default: `800`).

- height:

  The height of the video (default: `600`).

- framerate:

  The video framerate in frame/sec (default: `1`).

- res:

  The video resolution (default: `150`).

- quiet:

  A Boolean flag to enable/disable the messages.

## Value

The name of the produced video file path.

## Details

This function builds a video of the snapshots. It is available only if
the package `av` is installed.

## See also

`vignette("video")`,
[`plot_tissue()`](https://caravagnalab.github.io/ProCESS/1.3/reference/plot_tissue.md)
