"""
    AffineDomain2D(grid::CartesianDomain2D, map; units = units(grid))

A grid rotated, reflected, sheared, scaled or shifted with respect to the physical frame: the
point of the sample `[j, i]` is `map(grid[j, i])`. `grid` describes the samples in their own
frame (for example the pixels of a camera, in mm), and `map` takes a point of this frame to the
physical frame. `map` is an `AffineMap`, `LinearMap` or `Translation` of
CoordinateTransformations.jl; it is stored as an `AffineMap` with static 2×2 linear part, so
that indexing allocates nothing.

Arrays on the domain have the size and index order of `grid`.

`units` are those of the physical frame, one unit for both axes (a rotation mixes them); by
default the units of `grid`, which is right when `map` has no dimension (a rotation, a
reflection, a stretch).

```julia
using SampledDomains, CoordinateTransformations, StaticArrays

grid = make_centered_domain2D(640, 480, 5.5e-3)          # camera pixels, mm
θ = deg2rad(30)
rot = LinearMap(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])  # sensor rotated by 30°
dom = AffineDomain2D(grid, Translation(0.1, 0.0) ∘ rot)     # and shifted by 0.1 mm in x
dom[1, 1] == rot(grid[1, 1]) + SVector(0.1, 0.0)
```
"""
struct AffineDomain2D{T,D<:CartesianDomain2D,M<:AffineMap,U<:NamedTuple{(:x, :y)}} <:
       AbstractDomain{SVector{2,T},2}
    grid::D
    map::M
    units::U
    function AffineDomain2D(
        grid::CartesianDomain2D, map::Union{AffineMap,LinearMap,Translation};
        units=grid.units,
    )
        m = _affine2(map)
        T = promote_type(eltype(eltype(grid)), eltype(m.linear), eltype(m.translation))
        u = _units(units)
        u.x == u.y || throw(ArgumentError(
            "an affine domain has one unit for both axes, got $(u.x) and $(u.y)"))
        return new{T,typeof(grid),typeof(m),typeof(u)}(grid, m, u)
    end
end

_smatrix2(a::Number) = SMatrix{2,2}(a, zero(a), zero(a), a)
_smatrix2(a::AbstractMatrix) = SMatrix{2,2}(a)

_affine2(m::AffineMap) = AffineMap(_smatrix2(m.linear), SVector{2}(m.translation))
function _affine2(m::LinearMap)
    l = _smatrix2(m.linear)
    return AffineMap(l, zero(SVector{2,eltype(l)}))
end
function _affine2(m::Translation)
    t = SVector{2}(m.translation)
    return AffineMap(_smatrix2(one(eltype(t))), t)
end

Base.size(dom::AffineDomain2D) = size(dom.grid)

Base.@propagate_inbounds function Base.getindex(dom::AffineDomain2D{T}, j::Int, i::Int) where {T}
    @boundscheck checkbounds(dom, j, i)
    return SVector{2,T}(dom.map(@inbounds dom.grid[j, i]))
end

Base.:(==)(a::AffineDomain2D, b::AffineDomain2D) =
    a.grid == b.grid && a.map.linear == b.map.linear &&
    a.map.translation == b.map.translation && a.units == b.units

units(dom::AffineDomain2D) = dom.units
setunits(dom::AffineDomain2D, u) = AffineDomain2D(dom.grid, dom.map; units=u)

function convertunits(dom::AffineDomain2D, u)
    new = _units(u)
    f = factor(dom.units.x, new.x)
    return AffineDomain2D(dom.grid, AffineMap(dom.map.linear * f, dom.map.translation * f); units=new)
end

pointof(dom::AffineDomain2D, j::Real, i::Real) = dom.map(pointof(dom.grid, j, i))
indexof(dom::AffineDomain2D, p) = indexof(dom.grid, inv(dom.map)(SVector(p[1], p[2])))

"""
    dom * s
    s * dom

The affine domain with the physical coordinates multiplied by the number `s`; the grid is
unchanged. The units are dropped, as for a [`CartesianDomain2D`](@ref).
"""
Base.:*(dom::AffineDomain2D, s::Real) =
    AffineDomain2D(dom.grid, AffineMap(dom.map.linear * s, dom.map.translation * s); units=nothing)
Base.:*(s::Real, dom::AffineDomain2D) = dom * s
Base.:/(dom::AffineDomain2D, s::Real) = dom * inv(s)

function Base.show(io::IO, dom::AffineDomain2D)
    return print(io, "AffineDomain2D(", dom.grid, ", ", dom.map, ")")
end

function Base.show(io::IO, ::MIME"text/plain", dom::AffineDomain2D{T}) where {T}
    ny, nx = size(dom)
    println(io, ny, "×", nx, " AffineDomain2D{", T, "} with points map(grid[j, i]):")
    println(io, "  grid: ", dom.grid)
    return print(io, "  map:  ", dom.map)
end
