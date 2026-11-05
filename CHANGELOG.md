# Changelog

All notable changes to SampledDomains.jl are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the package follows
[semantic versioning](https://semver.org/).

## [0.2.0] — unreleased

### Changed (breaking)
- A domain is a lazy array of the points of its samples: `AbstractDomain{T,N} <: AbstractArray{T,N}`,
  and `CartesianDomain2D <: AbstractMatrix{SVector{2,T}}`. `dom[j, i]`, `dom[k]` and the
  iteration return the point as `SVector(x, y)` (with `p.x`, `p.y`), instead of the vector
  `[y, x]`. Indexing allocates nothing. `map(p -> f(p.x, p.y), dom)` gives the array of the
  values of `f` on the domain.
- `CartesianDomain2D{T,TX,TY}` and `SampledDomain{T,N,A,D}` are parametric, with concretely typed
  fields; code reading `dom.xrange`, `dom.yrange` is now type-stable.
- `size` and `length` no longer find the ranges by reflection; each domain type defines `size`.
- `getranges(dom)` returns a tuple `(xrange, yrange)` instead of a vector.
- `SampledDomain(vals, dom)` with sizes that differ throws a `DimensionMismatch` (before:
  `ErrorException`), and accepts any `AbstractArray` of values.
- `dualRange` is `fftshift(fftfreq(n * q, q / step))` of AbstractFFTs. Each frequency is now
  correct to one unit in the last place; before, the error was about one unit in the last place
  of the largest frequency, i.e. up to tens of units in the last place for small frequencies.
  Results that use the dual grid change in the last bits.
- Requires Julia 1.10.

### Added
- `AffineDomain2D(grid, map)`: a grid rotated, reflected, scaled or shifted with respect to the
  physical frame, with the map an `AffineMap`, `LinearMap` or `Translation` of
  CoordinateTransformations.jl.
- `dualDomain(dom, q::Int)` for the same upsampling in both directions.
- `show` methods: a domain prints as its ranges, not as a matrix of points.
- Documentation: the motivation of the convention, arrays of more than two dimensions,
  rotated grids.
- Tests: centred domains, scaling, `dualRange` for even and odd lengths and its accuracy
  (one unit in the last place, exact zero), Aqua.jl, allocation and inference of indexing.
- Dependencies: StaticArrays, CoordinateTransformations, AbstractFFTs.

## [0.1.0]

First version: `CartesianDomain2D`, `make_centered_domain2D`, `dualRange`, `dualDomain`,
`SampledDomain`.
