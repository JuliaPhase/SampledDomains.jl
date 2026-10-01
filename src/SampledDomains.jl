"""
The goal of the module is to attach physical coordinates to array elements, so an array entry could be interpreted as a sample of some function ``f(x, y)``.

Arrays follow the Julia image convention: the **first** index runs along `y` and the **second** along `x`, i.e.
`a[j, i]` is ``f(x_i, y_j)`` and `size(dom) == (length(yrange), length(xrange))`.
Row 1 corresponds to `yrange[1]`. To plot such an array with Makie's `heatmap`
(first index horizontal), transpose it: `heatmap(xrange, yrange, a')`.
"""
module SampledDomains

abstract type AbstractDomain end

import Base.:*
import Base: size, length, axes, getindex, ndims
import Base.iterate, Base.IteratorSize, Base.IndexStyle

export CartesianDomain2D, make_centered_domain2D, dualRange, dualDomain

getranges(dom::AbstractDomain) = [getfield(dom, n) for n in fieldnames(typeof(dom))]

length(dom::AbstractDomain) = prod(length.(getranges(dom)))
size(dom::AbstractDomain) = tuple(length.(reverse(getranges(dom)))...) # we need reverse to follow column major rule
# size(dom::AbstractDomain, d) = length(getfield(dom,d))
size(dom::AbstractDomain, d) = size(dom)[d]
ndims(dom::AbstractDomain) = length(fieldnames(typeof(dom)))

# TODO rewrite as parametric type
"""
    CartesianDomain2D(xrange, yrange)

Rectangular grid with sample coordinates `xrange` and `yrange`.
Arrays defined on it have size `(length(yrange), length(xrange))`,
so that `a[j, i]` is the value at `(xrange[i], yrange[j])`.

```julia
dom = CartesianDomain2D(-1:0.5:1, -1:0.5:1)
a = [x + 2y for y in dom.yrange, x in dom.xrange]   # a[j, i] = f(x[i], y[j])
```
"""
struct CartesianDomain2D <: AbstractDomain
    xrange::AbstractRange
    yrange::AbstractRange
end

*(dom::CartesianDomain2D, s::Real) = CartesianDomain2D(dom.xrange .* s, dom.yrange .* s)
*(s::Real, dom::CartesianDomain2D) = dom * s

"""
    dom[j, i]
    dom[k]

Coordinates of the sample with array index `(j, i)` (or linear index `k`), returned
as `[y, x]` — the same order as when iterating over `dom`.
"""
function getindex(dom::CartesianDomain2D, I::Vararg{Int,2})
    return collect(x[i] for (x, i) in zip(reverse(getranges(dom)), I))
end

## Linear indexing follows Julia's column-major order for an array of size `size(dom)`,
## i.e. it is consistent with `dom[j, i]` and with iteration over `dom`.
function Base.getindex(dom::CartesianDomain2D, i::Int)
    1 <= i <= length(dom) || throw(BoundsError(dom, i))
    xc, yc = divrem(i - 1, length(dom.yrange))
    return [dom.yrange[yc + 1], dom.xrange[xc + 1]]
end

Base.IndexStyle(::Type{<:CartesianDomain2D}) = IndexCartesian()

function iterate(dom::CartesianDomain2D, state)
    return iterate(Iterators.product(reverse(getranges(dom))...), state)
end
iterate(dom::CartesianDomain2D) = iterate(Iterators.product(reverse(getranges(dom))...))
# axes(dom::CartesianDomain2D) = axes(Iterators.product(getranges(dom)...))
# ndims(dom::CartesianDomain2D) = ndims(Iterators.product(getranges(dom)...))
IteratorSize(dom::CartesianDomain2D) = IteratorSize(Iterators.product(getranges(dom)...))

function make_centered_domain2D(xlength, ylength, pixelsizex, pixelsizey)
    xrange = ((1:xlength) .- (1 + xlength) / 2) .* pixelsizex
    yrange = ((1:ylength) .- (1 + ylength) / 2) .* pixelsizey
    return CartesianDomain2D(xrange, yrange)
end

make_centered_domain2D(xlength, ylength, pixelsize) =
    make_centered_domain2D(xlength, ylength, pixelsize, pixelsize)

"""
    dualRange(xrange::AbstractRange, q::Int =1)

Construct range in Fourier-transform-dual domain with upsampling factor `q`. If sampling in `x` has step size Δx, the dual
domain is sampled in interval (-q/2Δx, q/2Δx).
"""
function dualRange(xrange::AbstractRange, q::Int=1)
    len = length(xrange) * q
    st = step(xrange) / q
    return UnitRange(-floor(Int, len / 2), ceil(Int, len / 2) - 1) / (st * len)
end

function dualDomain(dom::CartesianDomain2D, q::Tuple{Int,Int}=(1, 1))
    kxrange = dualRange(dom.xrange, q[1])
    kyrange = dualRange(dom.yrange, q[2])
    return CartesianDomain2D(kxrange, kyrange)
end

"""
    SampledDomain(vals::Array, dom::AbstractDomain)
    SampledDomain(f::Function, dom::AbstractDomain)

Contains two fields: `vals` and `dom` representing sampled values of function `f` on `dom`ain.

When constructed from a function, `f` is called with the coordinates in the order of the
domain ranges, e.g. `f(x, y)` for a [`CartesianDomain2D`](@ref), and `vals` has the array layout
described there: `vals[j, i] = f(x[i], y[j])`.
"""
struct SampledDomain
    vals::Array
    dom::AbstractDomain
    function SampledDomain(vals::Array, dom::AbstractDomain)
        return if size(vals) != size(dom)
            throw(ErrorException("Different size of the domain and values array"))
        else
            new(vals, dom)
        end
    end
end

function SampledDomain(f::Function, dom::AbstractDomain)
    ## the product runs over the reversed ranges (first array index = last coordinate),
    ## so reverse each tuple back to call `f` in the natural order (x, y)
    return SampledDomain(
        map(t -> f(reverse(t)...), Iterators.product(reverse(getranges(dom))...)), dom
    )
end

end
