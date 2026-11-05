```@meta
CurrentModule = SampledDomains
```

# [Guide](@id guide)

## Grids and their points

A [`CartesianDomain2D`](@ref) is defined by two ranges. It is an `AbstractMatrix` of points
`SVector{2}`: indexing computes the point and allocates nothing.

```julia
using SampledDomains

dom = CartesianDomain2D(1:3, 10:10:20)   # 3 samples in x, 2 in y
size(dom)              # (2, 3)
dom[2, 3]              # [3, 20]: x = 3, y = 20
dom[2, 3].y            # 20
collect(dom)           # 2×3 Matrix of the points
dom.xrange, dom.yrange # the ranges
```

Everything that works on arrays works on a domain: iteration, `map`, broadcasting,
`CartesianIndices`, linear indexing (column-major, as for the arrays on it).

```julia
r2 = map(p -> p.x^2 + p.y^2, dom)                     # an Array of size(dom)
ap = map(p -> hypot(p.x, p.y) <= 15, dom)             # a mask
```

A domain prints as its ranges:

```julia
julia> dom
2×3 CartesianDomain2D{Int64} with points (x, y):
  x: 1:3
  y: 10:10:20
```

Multiplying by a number scales the coordinates, e.g. from pixels to millimetres:
`dom * 5.5e-3`.

## Centred grids

[`make_centered_domain2D`](@ref)`(nx, ny, pixelsize)` gives a grid of `nx × ny` samples with
the given step, centred on zero. For an even number of samples zero lies between the two
middle samples, for an odd number at the middle sample.

```julia
pupil = make_centered_domain2D(256, 256, 10e-3)   # 256 × 256 samples of 10 µm, in mm
```

## Fourier-dual grids

For a grid with step `Δx` and `n` samples, the discrete Fourier transform gives the spatial
frequencies with step `1 / (n Δx)`. [`dualRange`](@ref) and [`dualDomain`](@ref) return them in
the order of `fftshift`, with zero frequency at index `n ÷ 2 + 1`, and with an optional
upsampling factor `q` (zero padding of the field to `q n` samples):

```julia
fdom = dualDomain(pupil, 2)          # frequencies in 1/mm, upsampled twice
fdom * (λ * f)                       # coordinates in the focal plane of a lens of focal length f
```

The values agree with `fftshift(fftfreq(q * n, q / Δx))` of AbstractFFTs up to rounding.

## Rotated, reflected and shifted grids

When the grid of the samples is rotated, reflected or shifted with respect to the physical
frame (a camera with a rotated sensor, a beam on a mirror at oblique incidence), an
[`AffineDomain2D`](@ref) combines the grid in its own frame with a map to the physical frame,
from CoordinateTransformations.jl:

```julia
using CoordinateTransformations, StaticArrays

grid = make_centered_domain2D(640, 480, 5.5e-3)              # camera pixels, mm
θ = deg2rad(30)
rot = LinearMap(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])
dom = AffineDomain2D(grid, Translation(0.1, 0.0) ∘ rot)       # rotated by 30°, shifted
dom[1, 1] == rot(grid[1, 1]) + SVector(0.1, 0.0)
```

Arrays on an affine domain have the size and the index order of its grid; functions evaluated
on it see the physical coordinates. The map is stored as an `AffineMap` with a static matrix,
so its inverse `inv(dom.map)` takes physical points back to the frame of the grid.

## Values with their domain

[`SampledDomain`](@ref) keeps an array of values together with its domain and checks that
their sizes agree. Constructed from a function, it calls `f(x, y)` at every point:

```julia
s = SampledDomains.SampledDomain((x, y) -> x + 2y, dom)
s.vals[j, i] == dom.xrange[i] + 2dom.yrange[j]
```
