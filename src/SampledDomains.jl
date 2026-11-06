"""
SampledDomains attaches physical axes to arrays of data (camera images, pupils, PSFs), so that
it is always clear where each element lies.

A domain is a lazy array of the points of the samples: `dom[j, i]` is the point `(x, y)` of the
sample `a[j, i]` of an array `a` on the domain, returned as an `SVector` with the properties
`x` and `y`. Nothing is stored but the description of the grid and the units.

Arrays keep the index order of images, `a[j, i]` with row `j` and column `i`, and by default
a domain is oriented as an image read from a file: `x` to the right, `y` up, row 1 at the top
([`imagedomain`](@ref)). Then an image loaded from a file is at once the array of the values
on the domain, and the REPL prints an array as the picture looks.
"""
module SampledDomains

using StaticArrays: SVector, SMatrix
using CoordinateTransformations: AffineMap, LinearMap, Translation
using AbstractFFTs: fftfreq, fftshift

export AbstractDomain, CartesianDomain2D, AffineDomain2D
export imagedomain, make_centered_domain2D, dualRange, dualDomain
export LengthUnit, wavelengthunit, units, setunits, convertunits, inwavelengths, axislabels
export pixelsize, pointof, indexof
export imagej_calibration, imagej_domain, save_imagej, load_imagej

"""
    AbstractDomain{T,N} <: AbstractArray{T,N}

A domain: the array of the points of the samples of a grid, of the size of the arrays sampled
on it. A subtype defines

- `size(dom)`;
- `getindex(dom, I::Vararg{Int,N})`, the point of the sample `I`, as an `SVector`;
- `units(dom)`, the units of the coordinates;

and gets iteration, `map`, broadcasting and linear indexing from the `AbstractArray`
interface. Two-dimensional domains have points `SVector{2,T}` in the order `(x, y)`, and the
index order `(j, i)` of images: `dom[j, i]` is the point of `a[j, i]`.

Take the coordinates of a point by name, `p.x`, `p.y`, or `for (; x, y) in dom`, not by
position: such code does not depend on the order of the coordinates (in version 0.1 the point
was `[y, x]`).
"""
abstract type AbstractDomain{T,N} <: AbstractArray{T,N} end

include("units.jl")
include("cartesian.jl")
include("affine.jl")
include("sampled.jl")
include("imagej.jl")

end
