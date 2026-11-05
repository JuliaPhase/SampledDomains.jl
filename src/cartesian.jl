"""
    CartesianDomain2D(xrange, yrange)

Rectangular grid with sample coordinates `xrange` and `yrange`. Arrays defined on it have size
`(length(yrange), length(xrange))`, so that `a[j, i]` is the value at `(xrange[i], yrange[j])`.
The domain is the array of these points: `dom[j, i] == SVector(xrange[i], yrange[j])`, with
`dom[j, i].x == xrange[i]` and `dom[j, i].y == yrange[j]`.

```julia
dom = CartesianDomain2D(-1:0.5:1, -1:0.5:0)
a = map(p -> p.x + 2p.y, dom)            # a[j, i] = f(x[i], y[j]), size(a) == (3, 5)
a == [x + 2y for y in dom.yrange, x in dom.xrange]
```
"""
struct CartesianDomain2D{T,TX<:AbstractRange,TY<:AbstractRange} <:
       AbstractDomain{SVector{2,T},2}
    xrange::TX
    yrange::TY
    function CartesianDomain2D(xrange::AbstractRange, yrange::AbstractRange)
        T = promote_type(eltype(xrange), eltype(yrange))
        return new{T,typeof(xrange),typeof(yrange)}(xrange, yrange)
    end
end

Base.size(dom::CartesianDomain2D) = (length(dom.yrange), length(dom.xrange))

Base.@propagate_inbounds function Base.getindex(
    dom::CartesianDomain2D{T}, j::Int, i::Int
) where {T}
    @boundscheck checkbounds(dom, j, i)
    return @inbounds SVector{2,T}(dom.xrange[i], dom.yrange[j])
end

"""
    getranges(dom::CartesianDomain2D) -> (xrange, yrange)
"""
getranges(dom::CartesianDomain2D) = (dom.xrange, dom.yrange)

Base.:(==)(a::CartesianDomain2D, b::CartesianDomain2D) =
    a.xrange == b.xrange && a.yrange == b.yrange

"""
    dom * s
    s * dom

The domain with the coordinates multiplied by the number `s`.
"""
Base.:*(dom::CartesianDomain2D, s::Real) = CartesianDomain2D(dom.xrange .* s, dom.yrange .* s)
Base.:*(s::Real, dom::CartesianDomain2D) = dom * s

function Base.show(io::IO, dom::CartesianDomain2D)
    return print(io, "CartesianDomain2D(", dom.xrange, ", ", dom.yrange, ")")
end

function Base.show(io::IO, ::MIME"text/plain", dom::CartesianDomain2D{T}) where {T}
    ny, nx = size(dom)
    println(io, ny, "×", nx, " CartesianDomain2D{", T, "} with points (x, y):")
    println(io, "  x: ", dom.xrange)
    return print(io, "  y: ", dom.yrange)
end

"""
    make_centered_domain2D(xlength, ylength, pixelsize)
    make_centered_domain2D(xlength, ylength, pixelsizex, pixelsizey)

A [`CartesianDomain2D`](@ref) of `xlength × ylength` samples with the given steps, centred on
zero: for an even length zero lies between two samples, for an odd one at the middle sample.
"""
function make_centered_domain2D(xlength, ylength, pixelsizex, pixelsizey)
    xrange = ((1:xlength) .- (1 + xlength) / 2) .* pixelsizex
    yrange = ((1:ylength) .- (1 + ylength) / 2) .* pixelsizey
    return CartesianDomain2D(xrange, yrange)
end

make_centered_domain2D(xlength, ylength, pixelsize) =
    make_centered_domain2D(xlength, ylength, pixelsize, pixelsize)

"""
    dualRange(xrange::AbstractRange, q::Int = 1)

Range in the Fourier-dual domain with upsampling factor `q`. If the sampling in `x` has step
Δx, the dual range covers the interval `[-q/2Δx, q/2Δx)` with `length(xrange) * q` samples
and step `1 / (length(xrange) Δx)`, in the order of `fftshift`: zero frequency at index
`length ÷ 2 + 1`.

It is `fftshift(fftfreq(length(xrange) * q, q / step(xrange)))` of AbstractFFTs: each
frequency is computed as `k * (fs / n)`, correctly rounded to within one unit in the last
place, and zero is exact.
"""
function dualRange(xrange::AbstractRange, q::Int=1)
    return fftshift(fftfreq(length(xrange) * q, q / step(xrange)))
end

"""
    dualDomain(dom::CartesianDomain2D, q = (1, 1))

The Fourier-dual [`CartesianDomain2D`](@ref) of `dom`, with upsampling factors `q = (qx, qy)`
(or one integer for both), built from [`dualRange`](@ref).
"""
function dualDomain(dom::CartesianDomain2D, q::Tuple{Int,Int}=(1, 1))
    kxrange = dualRange(dom.xrange, q[1])
    kyrange = dualRange(dom.yrange, q[2])
    return CartesianDomain2D(kxrange, kyrange)
end

dualDomain(dom::CartesianDomain2D, q::Int) = dualDomain(dom, (q, q))
