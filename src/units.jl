"""
    LengthUnit(name::Symbol, metres, power = 1)

A unit of the coordinates of a domain: `name` for the labels, `metres` the length of one unit,
and `power` 1 for lengths or -1 for spatial frequencies (the dual domain). The units are a
label of the domain; the coordinates stay plain numbers.

Usually made from a symbol, `LengthUnit(:um)`, or by [`wavelengthunit`](@ref). Known symbols:
`:nm`, `:um` (also `:μm`, `:µm`), `:mm`, `:cm`, `:m`.
"""
struct LengthUnit
    name::Symbol
    metres::Float64
    power::Int
end

const _LENGTHS = (nm=1e-9, um=1e-6, mm=1e-3, cm=1e-2, m=1.0)
const _ALIASES = Dict(:μm => :um, :µm => :um, :micron => :um, :microns => :um)

LengthUnit(name::Symbol, metres::Real) = LengthUnit(name, metres, 1)

function LengthUnit(name::Symbol)
    n = get(_ALIASES, name, name)
    haskey(_LENGTHS, n) || throw(
        ArgumentError("unknown unit :$name; known: $(join(keys(_LENGTHS), ", ")), or wavelengthunit(λ, unit)"),
    )
    return LengthUnit(n, _LENGTHS[n], 1)
end

"""
    wavelengthunit(λ, unit)

The unit "one wavelength" for the wavelength `λ` given in `unit` (a symbol or a
[`LengthUnit`](@ref)): coordinates in this unit are `x/λ`, dimensionless numbers. The value
of λ is kept, for the labels and for the conversion back.

```julia
dom = imagedomain((64, 64); pixelsize=0.5, units=:um)
convertunits(dom, wavelengthunit(0.633, :um))     # coordinates x/λ
```
"""
wavelengthunit(λ::Real, unit) = LengthUnit(:λ, λ * _unit(unit).metres, 1)

_unit(u::LengthUnit) = u
_unit(u::Symbol) = LengthUnit(u)
_unit(::Nothing) = nothing

"""
    inv(u::LengthUnit)

The unit of the spatial frequencies dual to `u` (cycles per unit).
"""
Base.inv(u::LengthUnit) = LengthUnit(u.name, u.metres, -u.power)

Base.:(==)(a::LengthUnit, b::LengthUnit) =
    a.name == b.name && a.metres == b.metres && a.power == b.power

const _SUPERSCRIPT_MINUS1 = "⁻¹"

"""
    unitlabel(u)

The unit as text for an axis label: `"μm"`, `"λ"`, `"μm⁻¹"`; `""` for `nothing` (pixels or
dimensionless).
"""
unitlabel(::Nothing) = ""
function unitlabel(u::LengthUnit)
    s = u.name === :um ? "μm" : string(u.name)
    return u.power == 1 ? s : s * _SUPERSCRIPT_MINUS1
end

function Base.show(io::IO, u::LengthUnit)
    print(io, unitlabel(u))
    if u.name === :λ
        print(io, " (λ = ", _pretty(u.metres), ")")
    end
end

function _pretty(metres)
    for n in (:nm, :um, :mm, :m)
        v = metres / _LENGTHS[n]
        1 <= abs(v) < 1000 && return string(round(v; sigdigits=6), " ", n === :um ? "μm" : n)
    end
    return string(metres, " m")
end

# the units of a domain: a NamedTuple (x = ux, y = uy) of LengthUnit or nothing
_units(u::NamedTuple{(:x, :y)}) = (x=_unit(u.x), y=_unit(u.y))
_units(u::Union{Symbol,LengthUnit,Nothing}) = (x=_unit(u), y=_unit(u))

"""
    factor(from::LengthUnit, to::LengthUnit)

The number by which coordinates in `from` are multiplied to be in `to`. The units must have
the same power (both lengths or both frequencies).
"""
function factor(from::LengthUnit, to::LengthUnit)
    from.power == to.power ||
        throw(ArgumentError("cannot convert $(unitlabel(from)) to $(unitlabel(to))"))
    return from.power == 1 ? from.metres / to.metres : to.metres / from.metres
end
factor(::Nothing, ::Nothing) = 1
factor(::Nothing, to::LengthUnit) =
    throw(ArgumentError("the domain has no units; declare them first with setunits"))
factor(::LengthUnit, ::Nothing) =
    throw(ArgumentError("cannot convert to no units; use setunits(dom, nothing) to drop them"))

const _EXP10 = (nm=-9, um=-6, mm=-3, cm=-2, m=0)

# coordinates `r` in `from` converted to `to`. Between decimal units the range is multiplied
# or divided by an exact power of ten, so that 19.25 μm becomes 0.01925 mm, not
# 0.019249999999999996 mm.
function _rescale(r, from::LengthUnit, to::LengthUnit)
    if haskey(_EXP10, from.name) && haskey(_EXP10, to.name)
        from.power == to.power ||
            throw(ArgumentError("cannot convert $(unitlabel(from)) to $(unitlabel(to))"))
        e = (_EXP10[from.name] - _EXP10[to.name]) * from.power
        return e >= 0 ? r .* exp10(e) : r ./ exp10(-e)
    end
    factor(from, to)                                      # checks the powers
    # x/λ: divide by λ in the units of x (0.633), rather than multiply by its inverse
    lengths = from.power == 1
    if to.name === :λ
        return lengths ? r ./ (to.metres / from.metres) : r .* (to.metres / from.metres)
    else
        return lengths ? r .* (from.metres / to.metres) : r ./ (from.metres / to.metres)
    end
end
_rescale(r, from, to) = r .* factor(from, to)
