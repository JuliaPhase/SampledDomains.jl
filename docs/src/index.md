```@meta
CurrentModule = SampledDomains
```

# SampledDomains

SampledDomains attaches physical coordinates to array elements: a *domain* says at which point
`(x, y)` each element of an array was sampled. It is the bottom layer of the
[Phase.jl](https://github.com/JuliaPhase/Phase.jl) packages, which use it for pupils, bases,
focal-plane grids and camera images.

!!! warning "Draft"
    This documentation describes version 0.2, which is in development. Version 0.1 returned
    the coordinates of a sample as `[y, x]`; see the [changelog](https://github.com/JuliaPhase/SampledDomains.jl/blob/main/CHANGELOG.md).

## Why a domain

Julia has two conventions for two-dimensional arrays, and they disagree:

- **matrices and images** (JuliaImages, arrays read with FileIO): the first index is the row,
  that is the vertical direction;
- **plots** (Makie): in `heatmap(x, y, A)` and `image(A)` the first index runs along the
  horizontal axis.

Code that mixes images, FFTs and plots therefore keeps transposing arrays, and the sign or the
orientation of a tilt, a decentre or a rotation is easily lost. SampledDomains fixes one
convention, the one of images, and makes it explicit: the domain travels with the arrays and
knows which index is `x` and which is `y`. The arrays themselves stay plain arrays, so that
they go into FFTs, linear algebra and solvers without wrappers. See
[Conventions](@ref conventions).

## Quick start

```julia
using SampledDomains

dom = CartesianDomain2D(-1:0.1:1, -0.5:0.1:0.5)   # 21 samples in x, 11 in y
size(dom)                                         # (11, 21): rows are y, columns are x

p = dom[1, 2]                                     # the point of a[1, 2]
p.x, p.y                                          # (-0.9, -0.5)

a = map(p -> p.x + 2p.y, dom)                     # values of f(x, y) = x + 2y on the grid
a[1, 2] == -0.9 + 2 * (-0.5)
```

A domain is a lazy array of points, of the size of the arrays sampled on it: iteration,
`map`, broadcasting and indexing work as for any array, and nothing is stored but the ranges.

| Need | Use |
|---|---|
| a grid from coordinate ranges | [`CartesianDomain2D`](@ref) |
| a grid of pixels centred on zero | [`make_centered_domain2D`](@ref) |
| the Fourier-dual grid (spatial frequencies) | [`dualDomain`](@ref), [`dualRange`](@ref) |
| a rotated, reflected or shifted grid | [`AffineDomain2D`](@ref) |
| values of a function stored with their domain | [`SampledDomain`](@ref) |

The [Guide](@ref guide) shows each of them.
