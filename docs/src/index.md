```@meta
CurrentModule = SampledDomains
```

# SampledDomains

SampledDomains attaches physical axes to arrays of data — camera images, pupils, PSFs — so that
it is always clear where each element lies.

- **Oriented as an image by default.** An image imported with the standard libraries
  (`load` of FileIO) is shown with `x` to the right and `y` up, row 1 at the top. A domain made
  by [`imagedomain`](@ref) has the same orientation, so the loaded image is at once the array
  of the values on the domain, and the REPL prints it as the picture looks.
- **Origin chosen by the user:** the centre of the image (geometric, or the pixel of the zero
  frequency of an FFT) or the lower-left corner (of the image, or the centre of the corner
  pixel).
- **Units:** none by default (pixels). The pixel sizes Δx, Δy and their unit can be given at
  import or later; coordinates are converted to other units and to wavelengths, `x/λ`.
- **Physical coordinates for computations:** a function evaluated on the domain sees `x` and
  `y` in the units of the domain (a circle of radius 50 μm on an image); fractional indices
  are converted to points and back (distances between features, rasterisation).
- **Display and saving keep the axes:** Makie plots with the axes of the domain; TIFF files
  with the calibration of ImageJ (pixel size, unit, origin, direction of `y`).
- **Custom axes:** grids rotated, mirrored or stretched with respect to the physical frame
  ([`AffineDomain2D`](@ref)).

The arrays stay plain arrays: they go into FFTs, linear algebra and solvers without wrappers,
and the domain travels next to them.

!!! note "Version 0.2"
    Up to version 0.1, [`make_centered_domain2D`](@ref) gave `y` growing with the row number,
    so that a function evaluated on the domain was the vertical mirror of an image from a file,
    and the point of a sample was returned as `[y, x]`. Since 0.2, `y` points up and the point
    is `SVector(x, y)`. See the
    [changelog](https://github.com/JuliaPhase/SampledDomains.jl/blob/main/CHANGELOG.md).

## Quick start

An image of 480 rows and 640 columns from a camera with pixels of 5.5 μm, the origin at the
centre:

```@example quickstart
using SampledDomains
img = rand(480, 640)          # in practice: img = load("camera.tif")
dom = imagedomain(img; pixelsize=5.5, units=:um)
```

`dom[j, i]` is the point of the pixel `img[j, i]`, with coordinates `p.x` and `p.y`. The upper
left pixel has negative `x` and positive `y`:

```@example quickstart
p = dom[1, 1]
p.x, p.y
```

`map` evaluates a function at every point and returns an ordinary array of the size of the
image, here a disc of radius 50 μm:

```@example quickstart
disc = map(p -> hypot(p.x, p.y) <= 50, dom)
count(disc) * prod(pixelsize(dom)), π * 50^2      # its area in μm², and the exact one
```

## Contents

- [Images and axes](examples/ImagesAndAxes.md): an image from a file, its domain, how it is
  printed and plotted, stacks of images.
- [Origins and units](examples/OriginsAndUnits.md): the four origins, units and their
  conversion, a circle of 50 μm on an image, distances between features.
- [PSF in wavelengths](examples/PSF.md): the Fourier-dual domain, the focal plane in
  wavelengths, the direction in which a tilt moves the PSF.
- [ImageJ](examples/ImageJ.md): TIFF files with the calibration of ImageJ and Fiji.
- [Grids](examples/Grids.md): grids from ranges, rotated grids, values with their domain.
- [API](api.md): all types and functions.
