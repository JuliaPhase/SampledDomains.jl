"""
    CartesianDomain2D(xrange, yrange; units = nothing)

Rectangular grid with sample coordinates `xrange` and `yrange`. Arrays defined on it have size
`(length(yrange), length(xrange))`, so that `a[j, i]` is the value at `(xrange[i], yrange[j])`.
The domain is the array of these points: `dom[j, i] == SVector(xrange[i], yrange[j])`, with
`dom[j, i].x == xrange[i]` and `dom[j, i].y == yrange[j]`.

The ranges are taken as given. A domain oriented as an image, with `y` up, has a *descending*
`yrange`; [`imagedomain`](@ref) and [`make_centered_domain2D`](@ref) make such domains.

`units` labels the coordinates: `nothing` (pixels, or dimensionless), a symbol such as `:um`, a
[`LengthUnit`](@ref), or a NamedTuple `(x = …, y = …)` with a unit for each axis.

```julia
dom = CartesianDomain2D(-1:0.5:1, 0:-0.5:-1; units=:mm)
a = map(p -> p.x + 2p.y, dom)            # a[j, i] = f(x[i], y[j]), size(a) == (3, 5)
a == [x + 2y for y in dom.yrange, x in dom.xrange]
```
"""
struct CartesianDomain2D{T,TX<:AbstractRange,TY<:AbstractRange,U<:NamedTuple{(:x, :y)}} <:
       AbstractDomain{SVector{2,T},2}
    xrange::TX
    yrange::TY
    units::U
    function CartesianDomain2D(xrange::AbstractRange, yrange::AbstractRange; units=nothing)
        T = promote_type(eltype(xrange), eltype(yrange))
        u = _units(units)
        return new{T,typeof(xrange),typeof(yrange),typeof(u)}(xrange, yrange, u)
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
    a.xrange == b.xrange && a.yrange == b.yrange && a.units == b.units

# --- construction: image orientation and origin -------------------------------------------

"""
    imagedomain(sz::Tuple{Int,Int}; pixelsize = 1.0, origin = :center, units = nothing)
    imagedomain(a::AbstractMatrix; kwargs...)

A [`CartesianDomain2D`](@ref) for an array of size `sz = (rows, columns)` (or of the size of
`a`), oriented as an image read from a file: `x` grows to the right with the column `i`, `y`
grows *up*, i.e. decreases with the row `j`; row 1 is at the top.

- `pixelsize`: the steps `Δx`, `Δy` (positive), one number for both or a tuple `(Δx, Δy)`.
- `units`: the unit of the coordinates (see [`CartesianDomain2D`](@ref)).
- `origin`: where `(0, 0)` lies:
  - `:center`: the geometric centre of the array; for an even size between two pixels;
  - `:fftcenter`: the centre of the pixel `(rows ÷ 2 + 1, columns ÷ 2 + 1)`, where the
    zero frequency of `fftshift` lies;
  - `:lowerleft`: the lower-left corner of the image (the outer corner of the pixel);
  - `:lowerleftpixel`: the centre of the lower-left pixel;
  - a tuple `(j0, i0)`: a fractional index; pixel `(j, i)` has its centre at `(j, i)`, its
    corners at `j ± 0.5`, `i ± 0.5`.

```julia
img = load("camera.tif")                         # FileIO
dom = imagedomain(img; pixelsize=5.5, units=:um)
dom[1, 1]                                        # the upper-left pixel: x < 0, y > 0
```
"""
function imagedomain(
    sz::Tuple{Int,Int}; pixelsize=1.0, origin=:center, units=nothing
)
    m, n = sz
    Δx, Δy = _pair(pixelsize)
    (Δx > 0 && Δy > 0) || throw(ArgumentError("pixelsize must be positive, got $pixelsize"))
    j0, i0 = _originindex(origin, m, n)
    xrange = ((1:n) .- i0) .* Δx
    yrange = (j0 .- (1:m)) .* Δy
    return CartesianDomain2D(xrange, yrange; units)
end

imagedomain(a::AbstractMatrix; kwargs...) = imagedomain(size(a); kwargs...)

_pair(s::Real) = (float(s), float(s))
_pair(s::Tuple{Real,Real}) = (float(s[1]), float(s[2]))

function _originindex(origin::Symbol, m, n)
    origin === :center && return ((m + 1) / 2, (n + 1) / 2)
    origin === :fftcenter && return (float(m ÷ 2 + 1), float(n ÷ 2 + 1))
    origin === :lowerleft && return (m + 0.5, 0.5)
    origin === :lowerleftpixel && return (float(m), 1.0)
    throw(ArgumentError(
        "unknown origin :$origin; use :center, :fftcenter, :lowerleft, :lowerleftpixel or (j0, i0)",
    ))
end
_originindex(origin::Tuple{Real,Real}, m, n) = (float(origin[1]), float(origin[2]))

"""
    make_centered_domain2D(xlength, ylength, pixelsize; units = nothing)
    make_centered_domain2D(xlength, ylength, pixelsizex, pixelsizey; units = nothing)

A domain of `xlength × ylength` samples with the given steps, centred on zero and oriented as
an image (`y` up, descending `yrange`): `imagedomain((ylength, xlength); pixelsize, origin =
:center)`. For an even length zero lies between two samples, for an odd one at the middle
sample. Note the order of the arguments, `x` first.

Up to version 0.1 the `yrange` was ascending, so that `y` grew with the row.
"""
make_centered_domain2D(xlength, ylength, pixelsizex, pixelsizey; units=nothing) =
    imagedomain((ylength, xlength); pixelsize=(pixelsizex, pixelsizey), origin=:center, units)

make_centered_domain2D(xlength, ylength, pixelsize; units=nothing) =
    make_centered_domain2D(xlength, ylength, pixelsize, pixelsize; units)

# --- units ---------------------------------------------------------------------------------

"""
    units(dom) -> (x = ux, y = uy)

The units of the coordinates of a domain, each a [`LengthUnit`](@ref) or `nothing`.
"""
units(dom::CartesianDomain2D) = dom.units

"""
    setunits(dom, units)

The same domain, with the coordinates *declared* to be in `units` (the numbers are not
changed). `setunits(dom, nothing)` drops the units. To change the numbers, see
[`convertunits`](@ref).
"""
setunits(dom::CartesianDomain2D, u) = CartesianDomain2D(dom.xrange, dom.yrange; units=u)

"""
    convertunits(dom, units)

The same points in other units: a symbol (`:mm`), a [`LengthUnit`](@ref), e.g.
[`wavelengthunit`](@ref)`(0.633, :um)`, or a NamedTuple `(x = …, y = …)`. The domain must have
units.
"""
function convertunits(dom::CartesianDomain2D, u)
    new = _units(u)
    xrange = _rescale(dom.xrange, dom.units.x, new.x)
    yrange = _rescale(dom.yrange, dom.units.y, new.y)
    return CartesianDomain2D(xrange, yrange; units=new)
end

"""
    inwavelengths(dom, λ)

The domain with coordinates `x/λ`, `y/λ`, where `λ` is in the units of the domain:
`convertunits(dom, wavelengthunit(λ, units(dom).x))`. The axes of a PSF are then in wavelengths,
and the plots say so.
"""
function inwavelengths(dom::AbstractDomain, λ::Real)
    u = units(dom)
    u.x == u.y || throw(ArgumentError("the axes have different units, $(u.x) and $(u.y)"))
    u.x === nothing && throw(ArgumentError("the domain has no units; declare them first with setunits"))
    return convertunits(dom, wavelengthunit(λ, LengthUnit(u.x.name, u.x.metres, 1)))
end

"""
    axislabels(dom) -> (xlabel = "x, μm", ylabel = "y, μm")

Labels for the axes of a plot of an array on the domain, with the units; to be splatted into
a Makie `Axis`: `Axis(fig[1, 1]; axislabels(dom)...)`.
"""
function axislabels(dom::AbstractDomain)
    u = units(dom)
    lab(name, un) = un === nothing ? name : string(name, ", ", unitlabel(un))
    return (xlabel=lab("x", u.x), ylabel=lab("y", u.y))
end

"""
    pixelsize(dom::CartesianDomain2D) -> (Δx, Δy)

The steps of the grid, positive, in the units of the domain.
"""
pixelsize(dom::CartesianDomain2D) = (abs(step(dom.xrange)), abs(step(dom.yrange)))

# --- fractional indices ----------------------------------------------------------------------

"""
    pointof(dom, j, i)
    pointof(dom, (j, i))

The point of the (fractional) index `(j, i)`: the domain extended to real indices. For integer
indices it is `dom[j, i]` (up to rounding). Pixel `(j, i)` has its centre at the index
`(j, i)` and its corners at `j ± 0.5`, `i ± 0.5`. Used for the coordinates of features found in
an image (a centroid has fractional indices).
"""
function pointof(dom::CartesianDomain2D, j::Real, i::Real)
    x, y = dom.xrange, dom.yrange
    return SVector(first(x) + (i - 1) * step(x), first(y) + (j - 1) * step(y))
end
pointof(dom::AbstractDomain, ji::Union{Tuple{Real,Real},CartesianIndex{2}}) =
    pointof(dom, Tuple(ji)...)

"""
    indexof(dom, p) -> (j, i)

The fractional index of the point `p` (anything with `p[1] == x`, `p[2] == y`): the inverse of
[`pointof`](@ref). `round.(Int, indexof(dom, p))` is the pixel that contains `p`. Used for
rasterisation and to draw a point of known coordinates into an array.
"""
function indexof(dom::CartesianDomain2D, p)
    x, y = dom.xrange, dom.yrange
    return (1 + (p[2] - first(y)) / step(y), 1 + (p[1] - first(x)) / step(x))
end

# --- scaling and show -----------------------------------------------------------------------

"""
    dom * s
    s * dom
    dom / s

The domain with the coordinates multiplied (divided) by the number `s`. The units are dropped,
since a plain number does not tell what they become (dividing by a radius gives dimensionless
coordinates, multiplying a dual domain by `f λ` gives lengths); declare them again with
[`setunits`](@ref) if needed.
"""
Base.:*(dom::CartesianDomain2D, s::Real) = CartesianDomain2D(dom.xrange .* s, dom.yrange .* s)
Base.:*(s::Real, dom::CartesianDomain2D) = dom * s
Base.:/(dom::CartesianDomain2D, s::Real) = CartesianDomain2D(dom.xrange ./ s, dom.yrange ./ s)

_hasunits(u) = !(u.x === nothing && u.y === nothing)

function _showunits(io::IO, u)
    if u.x == u.y
        print(io, "units=", u.x === nothing ? "nothing" : repr(u.x))
    else
        print(io, "units=(x=", repr(u.x), ", y=", repr(u.y), ")")
    end
end

function Base.show(io::IO, dom::CartesianDomain2D)
    print(io, "CartesianDomain2D(", dom.xrange, ", ", dom.yrange)
    if _hasunits(dom.units)
        print(io, "; ")
        _showunits(io, dom.units)
    end
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", dom::CartesianDomain2D{T}) where {T}
    ny, nx = size(dom)
    println(io, ny, "×", nx, " CartesianDomain2D{", T, "} with points (x, y):")
    println(io, "  x: ", dom.xrange)
    print(io, "  y: ", dom.yrange)
    if _hasunits(dom.units)
        print(io, "\n  ")
        _showunits(io, dom.units)
    end
end

# --- Fourier-dual grid -----------------------------------------------------------------------

"""
    dualRange(xrange::AbstractRange, q::Int = 1)

Range in the Fourier-dual domain with upsampling factor `q`. If the sampling in `x` has step
Δx, the dual range covers the interval `[-q/2Δx, q/2Δx)` with `length(xrange) * q` samples
and step `1 / (length(xrange) Δx)`, in the order of `fftshift`: zero frequency at index
`length ÷ 2 + 1`.

A descending `xrange` (negative step, as `y` of an image) gives a descending dual range: the
frequency of the index `k` of an FFT is `k / (n Δx)`, with the sign of `Δx`.

It is `fftshift(fftfreq(length(xrange) * q, q / step(xrange)))` of AbstractFFTs: each
frequency is computed as `k * (fs / n)`, correctly rounded to within one unit in the last
place, and zero is exact.
"""
function dualRange(xrange::AbstractRange, q::Int=1)
    return fftshift(fftfreq(length(xrange) * q, q / step(xrange)))
end

"""
    dualDomain(dom::CartesianDomain2D, q = (1, 1))

The Fourier-dual [`CartesianDomain2D`](@ref) of `dom`, built from [`dualRange`](@ref), with
upsampling factors `q = (qx, qy)` (or one integer for both), and the inverse units (μm → μm⁻¹).
The upsampling keeps the frequency step `1 / (n Δx)` and makes the range `q` times wider: it is
the dual grid of an array over the same extent sampled `q` times finer. Zero padding, which
makes the frequency step finer, corresponds instead to the dual of a larger domain with the same
step.

For a domain oriented as an image (`y` up), the dual domain is oriented in the same way: the
`fftshift`ed FFT of an array on `dom` is an array on `dualDomain(dom)`, with positive
frequencies of `y` in the upper half.
"""
function dualDomain(dom::CartesianDomain2D, q::Tuple{Int,Int}=(1, 1))
    kxrange = dualRange(dom.xrange, q[1])
    kyrange = dualRange(dom.yrange, q[2])
    u = dom.units
    return CartesianDomain2D(kxrange, kyrange; units=(x=_inv(u.x), y=_inv(u.y)))
end

dualDomain(dom::CartesianDomain2D, q::Int) = dualDomain(dom, (q, q))

_inv(::Nothing) = nothing
_inv(u::LengthUnit) = inv(u)
