```@meta
CurrentModule = SampledDomains
```

# SampledDomains

SampledDomains attaches physical coordinates to array elements: a domain says at which point
`(x, y)` each element of an array was sampled.

## Why a domain

Julia has two conventions for two-dimensional arrays, and they disagree:

- **matrices and images** (JuliaImages, arrays read with FileIO): the first index is the row,
  that is the vertical direction;
- **plots** (Makie): in `heatmap(x, y, A)` and `image(A)` the first index runs along the
  horizontal axis.

Code that mixes images, FFTs and plots therefore keeps transposing arrays, and the sign or the
orientation of a tilt, a decentre or a rotation is easily lost. SampledDomains fixes one
convention, the one of images, and makes it explicit: the domain object travels with the
arrays and knows which index is `x` and which is `y`. The arrays themselves stay plain arrays,
so that they go into FFTs, linear algebra and solvers without wrappers.

## Array axis convention

Arrays defined on a [`CartesianDomain2D`](@ref) follow the image convention:
the **first** array index runs along `y` and the **second** along `x`.

```
a[j, i] = f(x[i], y[j])        size(a) == (length(yrange), length(xrange))
```

Row 1 corresponds to `yrange[1]`. This is the same index order as for images loaded
with `FileIO`/`Images` and the layout expected by FFT-based code, so arrays on a domain can be
combined with image data without transposition.

## A domain is an array of points

A domain is a lazy array of the points of the samples, of the same size as the arrays on it:
`dom[j, i]` is the point of `a[j, i]`, as an `SVector` in the natural order `(x, y)`, with the
properties `x` and `y`. Nothing is stored but the ranges, as for `CartesianIndices`.

```julia
using SampledDomains

dom = CartesianDomain2D(-1:0.1:1, -0.5:0.1:0.5)   # 21 samples in x, 11 in y
size(dom)                                         # (11, 21)
p = dom[1, 2]                                     # SVector(-0.9, -0.5)
p.x, p.y                                          # (-0.9, -0.5)

a = map(p -> p.x + 2p.y, dom)                     # a[j, i] = x[i] + 2y[j]
a == [x + 2y for y in dom.yrange, x in dom.xrange]
```

Iteration, `map`, broadcasting and linear indexing come from the `AbstractArray` interface.
`SampledDomain(f, dom)` stores the values of `f(x, y)` together with the domain:

```julia
s = SampledDomains.SampledDomain((x, y) -> x + 2y, dom)
s.vals[j, i] == dom.xrange[i] + 2dom.yrange[j]
```

Centred grids and the Fourier-dual grid (in the order of `fftshift`, with upsampling):

```julia
dom = make_centered_domain2D(256, 256, 10e-3)     # 256 × 256 pixels of 10 µm, centred on zero
fdom = dualDomain(dom, 2)                         # spatial frequencies, upsampled twice
```

## Rotated and reflected grids

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
dom[j, i] == rot(grid[j, i]) + SVector(0.1, 0.0)
```

Arrays on an affine domain have the size and the index order of its grid; functions evaluated
on it see the physical coordinates.

## More than two dimensions

The spatial axes come first, in the order `(y, x)`, followed by the other axes (frame, defocus,
mode number), as in JuliaImages. The domain describes the first two axes:

```julia
stack = zeros(size(dom)..., 5)          # five frames on dom
for frame in eachslice(stack; dims=3)   # each slice is an array on dom
    frame .= map(p -> p.x^2 + p.y^2, dom)
end
```

## Displaying arrays

Makie's `heatmap(A)` puts the **first** index on the horizontal axis, so it shows
`a` with `x` and `y` interchanged. To get `x` to the right and `y` up, transpose:

```julia
heatmap(dom.xrange, dom.yrange, a')    # axis values are shown, y increases upward
```

With `PhasePlots` this is done for you: `showarray(dom, a)` shows `a` with `x` to the right and `y` up,
real coordinates on the axes, and axes labelled `x` and `y`.

## Relation to images read from files

An image loaded from a file has row 1 at the **top**, so a Cartesian frame with `y` up
is `y_user = -ind1`. An array on a domain with an ascending `yrange` has row 1 at the
**bottom** of a y-up plot. The two therefore differ by a vertical flip. To work in the
image frame, either use a descending `yrange`, or flip the array (`reverse(a; dims=1)`)
before comparing. Both frames use the same index order, so no transposition is needed.

## API

```@index
```

```@autodocs
Modules = [SampledDomains]
```
