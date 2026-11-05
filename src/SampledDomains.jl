"""
The goal of the module is to attach physical coordinates to array elements, so an array entry
could be interpreted as a sample of some function ``f(x, y)``.

A domain is a lazy array of the points of the samples: `dom[j, i]` is the point `(x, y)` of the
sample `a[j, i]` of an array `a` on the domain, returned as an `SVector` with the properties
`x` and `y`. Nothing is stored but the description of the grid.

Arrays follow the Julia image convention: the **first** index runs along `y` and the
**second** along `x`, i.e. `a[j, i]` is ``f(x_i, y_j)`` and
`size(dom) == (length(yrange), length(xrange))`. Row 1 corresponds to `yrange[1]`. To plot
such an array with Makie's `heatmap` (first index horizontal), transpose it:
`heatmap(xrange, yrange, a')`.
"""
module SampledDomains

using StaticArrays: SVector, SMatrix
using CoordinateTransformations: AffineMap, LinearMap, Translation

export AbstractDomain, CartesianDomain2D, AffineDomain2D
export make_centered_domain2D, dualRange, dualDomain

"""
    AbstractDomain{T,N} <: AbstractArray{T,N}

A domain: the array of the points of the samples of a grid, of the size of the arrays sampled
on it. A subtype defines

- `size(dom)`;
- `getindex(dom, I::Vararg{Int,N})`, the point of the sample `I`, as an `SVector`;

and gets iteration, `map`, broadcasting and linear indexing from the `AbstractArray`
interface. Two-dimensional domains have points `SVector{2,T}` in the order `(x, y)`, and the
index order `(j, i)` of images: `dom[j, i]` is the point of `a[j, i]`.
"""
abstract type AbstractDomain{T,N} <: AbstractArray{T,N} end

include("cartesian.jl")
include("affine.jl")
include("sampled.jl")

end
