```@meta
CurrentModule = SampledDomains
```

# SampledDomains

SampledDomains attaches physical coordinates to the elements of an array: a domain says at
which point `(x, y)` each element was sampled. The Phase.jl packages use it for pupils, bases,
focal-plane grids and camera images.

!!! note "Draft"
    This documentation describes version 0.2, which is in development. In version 0.1 the
    point of a sample was returned as `[y, x]`; see the
    [changelog](https://github.com/JuliaPhase/SampledDomains.jl/blob/main/CHANGELOG.md).

## Why

Julia has two conventions for two-dimensional arrays. For matrices and images the first index
is the row, that is the vertical direction. For Makie plots the first index runs along the
horizontal axis. Code that mixes images, FFTs and plots keeps transposing arrays, and the sign
of a tilt or the direction of a rotation is easily lost.

SampledDomains takes the convention of images and makes it explicit: a domain travels with
the arrays and knows which index is `x` and which is `y`. The arrays stay plain arrays, so
they go into FFTs and solvers without wrappers. The page
[Arrays and pictures](examples/ArraysAndPictures.md) shows the difference in pictures.

## Quick start

A grid with 21 samples in `x` and 11 in `y`:

```@example quickstart
using SampledDomains
dom = CartesianDomain2D(-1:0.1:1, -0.5:0.1:0.5)
```

The domain has the size of the arrays on it: rows are `y`, columns are `x`.

```@example quickstart
size(dom)
```

`dom[j, i]` is the point of the sample `[j, i]`, with coordinates `p.x` and `p.y`:

```@example quickstart
p = dom[1, 2]
p.x, p.y
```

`map` evaluates a function at every point and returns an ordinary array:

```@example quickstart
a = map(p -> p.x + 2p.y, dom)
size(a)
```

## Contents

- [Arrays and pictures](examples/ArraysAndPictures.md): the index order, the direction of `y`,
  how to plot, stacks of arrays.
- [Grids](examples/Grids.md): grids from ranges, grids centred on zero, grids of spatial
  frequencies, rotated grids.
- [API](api.md): all types and functions.
