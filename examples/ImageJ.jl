# ```@meta
# CurrentModule = SampledDomains
# ```

using SampledDomains
using TiffImages
using CairoMakie
CairoMakie.activate!(; type="png")

# # ImageJ
#
# ImageJ and Fiji keep the calibration of an image in its TIFF file: the pixel size, the unit,
# the origin and the direction of `y`. SampledDomains writes and reads this calibration, so a
# domain survives a round trip through ImageJ, and an image calibrated in ImageJ comes with its
# domain. Writing needs TiffImages (`using TiffImages`); reading the calibration needs nothing.
#
# ## Saving
#
# An image of 60 × 80 pixels of 5.5 μm, origin at the lower-left corner:

dom = imagedomain((60, 80); pixelsize=5.5, origin=:lowerleft, units=:um)
img = map(p -> sinpi(p.x / 100) * cospi(p.y / 70), dom)
file = joinpath(mktempdir(), "calibrated.tif")
save_imagej(file, img, dom)

# The calibration as ImageJ will see it, from [`imagej_calibration`](@ref): the origin in
# pixels from the upper-left corner of the image, and `inverty`, since `y` of the domain points
# up while ImageJ counts pixel rows down:

imagej_calibration(dom)

# Opened in ImageJ (or Fiji), the image shows these coordinates in the status bar and in the
# measurements: the centre of the upper-left pixel at `(2.75, 327.25)` μm, as `dom[1, 1]`:

dom[1, 1]

# ## Loading
#
# [`load_imagej`](@ref) returns the image and its domain; [`imagej_domain`](@ref)`(file)`
# reads only the domain, for an image loaded otherwise, e.g. with `load` of FileIO:

img2, dom2 = load_imagej(file)
dom2 == dom

# The values come back as 32-bit floats, equal to the original ones up to the rounding to
# `Float32`:

maximum(abs.(Float64.(img2) .- img))

#-

@assert dom2 == dom && isapprox(Float64.(img2), img; rtol=1e-6)

# A file calibrated in ImageJ without "Invert Y coordinates" (in *Analyze ▸ Set Measurements*
# or *Image ▸ Properties*) has `y` pointing down, as ImageJ shows it then; its domain has an
# ascending `yrange`. A file without the calibration of ImageJ gives a domain in pixels,
# oriented as an image, with the origin at the centre.
#
# ## Details
#
# - ImageJ counts coordinates from the corner of a pixel: the pixel `(0, 0)` covers
#   `[0, 1] × [0, 1]`, its centre is `(0.5, 0.5)`. The pixel `[j, i]` of Julia has its centre at
#   `(i - 0.5, j - 0.5)` in ImageJ's pixel coordinates.
# - ImageJ stores the pixel size as a fraction, "pixels per unit"; it writes 5.5 μm as
#   181818/1000000 and reads back 5.5000055 μm. [`save_imagej`](@ref) writes the exact fraction
#   2/11.
# - With `inverty=true`, ImageJ takes `yorigin = 0` as "no origin given" and then puts `y = 0`
#   one pixel above the lower edge. A domain with `y` up and the origin on the upper edge of the
#   image would need exactly this value, so [`imagej_calibration`](@ref) writes
#   [`IMAGEJ_ZERO_YORIGIN`](@ref) `= 1e-9` pixels instead, invisible in ImageJ, and
#   [`imagej_domain`](@ref) reads it back as zero:

top = imagedomain((60, 80); pixelsize=5.5, origin=(0.5, 40.5), units=:um)
imagej_calibration(top).yorigin

# - Real values are saved as 32-bit floats, the type ImageJ works with. `save_imagej(file, a,
#   dom; float64=true)` saves 64-bit floats, read back exactly by [`load_imagej`](@ref); such
#   a file is for Julia: ImageJ 1.54 does not open it (only through *File ▸ Import ▸ Raw*, without
#   the calibration), and Fiji opens it through Bio-Formats as 32 bit, with the pixel size and
#   unit but without the origin and the direction of `y`. Units: nm, μm (written as
#   `micron`), mm, cm, m; a domain in wavelengths is written with the unit `lambda` and the
#   wavelength in metres. Different units of `x` and `y` are written as `unit` and `yunit`.
# - The translation between a domain and the calibration of ImageJ is in one file of the
#   package (`src/imagej.jl`); it was checked with ImageJ 1.54g.
