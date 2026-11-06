module SampledDomainsTiffImagesExt

using SampledDomains
using SampledDomains: imagej_description
using TiffImages
using TiffImages: DenseTaggedImage, ifds
using TiffImages: Gray, N0f8, N0f16      # ColorTypes, FixedPointNumbers through TiffImages

_gray(a::AbstractMatrix{<:Gray}) = a
_gray(a::AbstractMatrix{Float64}) = Gray{Float32}.(a)
_gray(a::AbstractMatrix{Float32}) = Gray{Float32}.(a)
_gray(a::AbstractMatrix{UInt8}) = Gray.(reinterpret(N0f8, a))
_gray(a::AbstractMatrix{UInt16}) = Gray.(reinterpret(N0f16, a))
_gray(a::AbstractMatrix{Bool}) = Gray{N0f8}.(a)
_gray(a::AbstractMatrix{<:Real}) = Gray{Float32}.(a)

# 1 / pixel size as a rational of UInt32, exact for sizes like 5.5 (2/11)
_resolution(p::Real) = rationalize(UInt32, 1 / p; tol=eps(1 / p))

function SampledDomains.save_imagej(path::AbstractString, a::AbstractMatrix,
    dom::CartesianDomain2D; float64::Bool=false)
    size(a) == size(dom) ||
        throw(DimensionMismatch("array of size $(size(a)) on a domain of size $(size(dom))"))
    cal = imagej_calibration(dom)
    g = float64 && eltype(a) <: Real ? Gray{Float64}.(a) : _gray(a)
    img = DenseTaggedImage(collect(g))
    ifd = ifds(img)
    ifd[TiffImages.IMAGEDESCRIPTION] = imagej_description(cal)
    ifd[TiffImages.XRESOLUTION] = _resolution(cal.pixelwidth)
    ifd[TiffImages.YRESOLUTION] = _resolution(cal.pixelheight)
    ifd[TiffImages.RESOLUTIONUNIT] = UInt16(1)          # none; the unit is in the description
    TiffImages.save(path, img)
    return path
end

function SampledDomains.load_imagej(path::AbstractString)
    return TiffImages.load(path), imagej_domain(path)
end

end
