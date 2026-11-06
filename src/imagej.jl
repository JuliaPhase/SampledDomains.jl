# The calibration of ImageJ (and Fiji): what a domain becomes in a TIFF file and back.
# Everything specific to ImageJ is here, so that a change of ImageJ is fixed in one place.
# Checked with ImageJ 1.54g (Phase/plots/SampledDomains/ImageJ_check).
#
# ImageJ counts coordinates from the corner of a pixel: pixel (0, 0) covers [0, 1]², its centre
# is (0.5, 0.5), so the pixel (j, i) of Julia (1-based) has its centre at (i - 0.5, j - 0.5).
# Scaled coordinates: x = (x_px - xorigin) * pw; y = (y_px - yorigin) * ph (y down), or with
# inverty=true y = (yorigin - y_px) * ph (y up). Quirk: with inverty=true and yorigin == 0
# ImageJ uses y = (H - 1 - y_px) * ph, as if the origin were not set.

"""
The `yorigin` written instead of zero for a domain with `y` up and the origin on the upper edge
of the image, since ImageJ takes `yorigin = 0` with `inverty=true` as "no origin" (it then puts
`y = 0` one pixel above the bottom edge). In pixels; invisible in ImageJ.
"""
const IMAGEJ_ZERO_YORIGIN = 1e-9

const _IMAGEJ_UNITNAMES = Dict(:nm => "nm", :um => "micron", :mm => "mm", :cm => "cm", :m => "m")
const _IMAGEJ_UNITS = Dict(
    "nm" => :nm, "nanometer" => :nm,
    "micron" => :um, "microns" => :um, "um" => :um, "µm" => :um, "μm" => :um,
    "mm" => :mm, "millimeter" => :mm, "cm" => :cm, "m" => :m, "meter" => :m,
)

"""
    imagej_calibration(dom::CartesianDomain2D) -> NamedTuple

The calibration of ImageJ for an image on `dom`: `pixelwidth`, `pixelheight` (positive),
`unit`, `yunit` (`nothing` if as `unit`), `xorigin`, `yorigin` (in pixels, from the corner of the
image), `inverty` (true for a domain with `y` up), `wavelength` (in metres, for a unit of
[`wavelengthunit`](@ref), else `nothing`).

The inverse is [`imagej_domain`](@ref). With `y` up and the origin on the upper edge,
`yorigin` is written as [`IMAGEJ_ZERO_YORIGIN`](@ref) instead of zero. ImageJ has no reversed
`x` axis, so a descending `xrange` is an error.
"""
function imagej_calibration(dom::CartesianDomain2D)
    x, y = dom.xrange, dom.yrange
    Δx, Δy = float(step(x)), float(step(y))
    Δx > 0 || throw(ArgumentError("ImageJ has no reversed x axis; the xrange descends"))
    inverty = Δy < 0
    ph = abs(Δy)
    xorigin = 0.5 - first(x) / Δx
    yorigin = inverty ? 0.5 + first(y) / ph : 0.5 - first(y) / ph
    if inverty && yorigin == 0
        yorigin = IMAGEJ_ZERO_YORIGIN
    end
    u = dom.units
    unit = _imagej_unitname(u.x)
    yunit = u.y == u.x ? nothing : _imagej_unitname(u.y)
    λ = _wavelength_metres(u.x, u.y)
    return (; pixelwidth=Δx, pixelheight=ph, unit, yunit, xorigin, yorigin, inverty,
        wavelength=λ)
end

_imagej_unitname(::Nothing) = "pixel"
function _imagej_unitname(u::LengthUnit)
    base = u.name === :λ ? "lambda" : get(_IMAGEJ_UNITNAMES, u.name, string(u.name))
    return u.power == 1 ? base : "1/" * base
end

function _wavelength_metres(ux, uy)
    for u in (ux, uy)
        u isa LengthUnit && u.name === :λ && return u.metres
    end
    return nothing
end

"""
    imagej_description(cal) -> String

The text of the ImageDescription tag of a TIFF file for the calibration `cal` of
[`imagej_calibration`](@ref), in the format that ImageJ writes and reads.
"""
function imagej_description(cal)
    io = IOBuffer()
    println(io, "ImageJ=1.54g")
    println(io, "unit=", cal.unit)
    cal.yunit === nothing || println(io, "yunit=", cal.yunit)
    println(io, "xorigin=", cal.xorigin)
    println(io, "yorigin=", cal.yorigin)
    cal.inverty && println(io, "inverty=true")
    cal.wavelength === nothing || println(io, "wavelength=", cal.wavelength)
    return String(take!(io))
end

"""
    parse_imagej_description(s) -> Dict{String,String}

The `key=value` lines of an ImageDescription tag written by ImageJ.
"""
function parse_imagej_description(s::AbstractString)
    d = Dict{String,String}()
    for line in split(s, '\n')
        k = findfirst('=', line)
        k === nothing && continue
        d[strip(line[1:prevind(line, k)])] = strip(line[nextind(line, k):end])
    end
    return d
end

"""
    imagej_domain(sz, description::AbstractString, pixelwidth, pixelheight) -> CartesianDomain2D

The domain of an image of size `sz = (rows, columns)` from the ImageDescription text of ImageJ
and the pixel sizes (`1 / XResolution`, `1 / YResolution`), as ImageJ itself interprets them,
including the case `yorigin = 0` with `inverty=true`. A file without `inverty=true` gives a
domain with `y` down (ascending `yrange`), as ImageJ shows it.

`imagej_domain(path)` reads the domain of a TIFF file.
"""
function imagej_domain(sz::Tuple{Int,Int}, description::AbstractString, pw::Real, ph::Real)
    m, n = sz
    d = parse_imagej_description(description)
    getnum(k) = haskey(d, k) ? parse(Float64, d[k]) : 0.0
    xorigin, yorigin = getnum("xorigin"), getnum("yorigin")
    inverty = get(d, "inverty", "false") == "true"
    λ = haskey(d, "wavelength") ? parse(Float64, d["wavelength"]) : nothing
    ux = _imagej_unit(get(d, "unit", "pixel"), λ)
    uy = haskey(d, "yunit") ? _imagej_unit(d["yunit"], λ) : ux
    xrange = ((1:n) .- (xorigin + 0.5)) .* float(pw)
    j0 = if !inverty
        yorigin + 0.5
    elseif yorigin == 0                      # the quirk of ImageJ
        m - 0.5
    elseif yorigin == IMAGEJ_ZERO_YORIGIN    # written by us for an origin on the upper edge
        0.5
    else
        yorigin + 0.5
    end
    yrange = inverty ? (j0 .- (1:m)) .* float(ph) : ((1:m) .- j0) .* float(ph)
    return CartesianDomain2D(xrange, yrange; units=(x=ux, y=uy))
end

function _imagej_unit(name::AbstractString, λ)
    s = strip(name)
    inverse = startswith(s, "1/")
    base = inverse ? s[3:end] : s
    u = if base in ("pixel", "pixels", "")
        nothing
    elseif base == "lambda"
        λ === nothing && throw(ArgumentError("unit lambda without the key wavelength="))
        LengthUnit(:λ, λ, 1)
    elseif haskey(_IMAGEJ_UNITS, base)
        LengthUnit(_IMAGEJ_UNITS[base])
    else
        @warn "unit \"$base\" of the TIFF file is not known; the domain has no units"
        nothing
    end
    return inverse && u !== nothing ? inv(u) : u
end

"""
    save_imagej(path, a::AbstractMatrix, dom::CartesianDomain2D)

Save `a` as a TIFF file with the calibration of `dom` in the format of ImageJ (pixel size, unit,
origin, direction of `y`; see [`imagej_calibration`](@ref)), readable by ImageJ and Fiji and by
[`load_imagej`](@ref). `Float64` values are saved as `Float32` (ImageJ reads 32-bit floats).
Needs TiffImages: `using TiffImages`.
"""
function save_imagej end

"""
    load_imagej(path) -> (img, dom)

Load a TIFF file and the domain of its ImageJ calibration (see [`imagej_domain`](@ref)). A file
without the calibration of ImageJ gives a domain with no units, oriented as an image, with the
origin at the centre. Needs TiffImages: `using TiffImages`.
"""
function load_imagej end

# --- the tags of a TIFF file, read directly ------------------------------------------------
# TiffImages 0.11 swaps numerator and denominator of RATIONAL tags in big-endian files (ImageJ
# writes big-endian files), so the calibration is read here, from the first IFD.

const _TIFF_WANTED = (256, 257, 270, 282, 283)   # width, height, description, x/y resolution

function _tiff_tags(path::AbstractString)
    b = read(path)
    length(b) >= 8 || throw(ArgumentError("$path is not a TIFF file"))
    le = b[1:2] == b"II"
    (le || b[1:2] == b"MM") || throw(ArgumentError("$path is not a TIFF file"))
    rd(T, o) = (v = reinterpret(T, b[o+1:o+sizeof(T)])[1]; le ? ltoh(v) : ntoh(v))
    rd(UInt16, 2) == 42 || throw(ArgumentError("$path is not a classic TIFF file"))
    off = Int(rd(UInt32, 4))
    tags = Dict{Int,Any}()
    for k in 0:Int(rd(UInt16, off))-1
        e = off + 2 + 12k
        tag, typ, cnt = Int(rd(UInt16, e)), Int(rd(UInt16, e + 2)), Int(rd(UInt32, e + 4))
        tag in _TIFF_WANTED || continue
        if typ == 2                                       # ASCII
            v = cnt <= 4 ? e + 8 : Int(rd(UInt32, e + 8))
            tags[tag] = String(b[v+1:v+cnt]) |> s -> rstrip(s, '\0')
        elseif typ == 3                                   # SHORT
            tags[tag] = Int(rd(UInt16, e + 8))
        elseif typ == 4                                   # LONG
            tags[tag] = Int(rd(UInt32, e + 8))
        elseif typ == 5                                   # RATIONAL
            v = Int(rd(UInt32, e + 8))
            tags[tag] = rd(UInt32, v) / rd(UInt32, v + 4)
        end
    end
    return tags
end

"""
    imagej_domain(path::AbstractString) -> CartesianDomain2D

The domain of the image in the TIFF file `path`, from its ImageJ calibration (see
[`imagej_domain`](@ref) with a description). A file without the calibration of ImageJ gives
`imagedomain(size)`: no units, oriented as an image, origin at the centre. The pixels are read
separately, e.g. with `load(path)` of FileIO; this needs no package.
"""
function imagej_domain(path::AbstractString)
    t = _tiff_tags(path)
    sz = (t[257], t[256])
    desc = get(t, 270, "")
    startswith(desc, "ImageJ=") || return imagedomain(sz)
    pw = 1 / get(t, 282, 1.0)
    ph = 1 / get(t, 283, 1.0)
    return imagej_domain(sz, desc, pw, ph)
end
