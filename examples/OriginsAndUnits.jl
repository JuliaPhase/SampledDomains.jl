# ```@meta
# CurrentModule = SampledDomains
# ```

using SampledDomains
using CairoMakie
using Random
CairoMakie.activate!(; type="png")

# # Origins and units
#
# ## Four origins
#
# [`imagedomain`](@ref) puts the origin `(0, 0)` at one of four places, chosen by the keyword
# `origin`. Here an image of 4 rows and 5 columns with pixels of size 1; each pixel is drawn as
# a square around its point, and the origin as a red dot:

origins = [:center, :fftcenter, :lowerleft, :lowerleftpixel]
fig = Figure(; size=(1000, 300))
for (k, o) in enumerate(origins)
    d = imagedomain((4, 5); origin=o)
    ax = Axis(fig[1, k]; title=":$o", aspect=DataAspect(), xlabel="x", ylabel="y",
        limits=(-3.2, 5.7, -2.7, 4.2))
    for p in d
        poly!(ax, Rect2f(p.x - 0.5, p.y - 0.5, 1, 1); color=:white, strokecolor=:gray40,
            strokewidth=1)
    end
    text!(ax, d[1, 1]...; text="[1,1]", align=(:center, :center), fontsize=10)
    scatter!(ax, [Point2f(0, 0)]; color=:red, markersize=10)
end
fig

# - `:center`: the geometric centre of the image; for an even number of pixels it lies between
#   two pixels (the default).
# - `:fftcenter`: the centre of the pixel `(rows ÷ 2 + 1, columns ÷ 2 + 1)`, where `fftshift`
#   puts the zero frequency. A pupil on such a domain goes into an FFT with `ifftshift`
#   without a phase ramp.
# - `:lowerleft`: the lower-left corner of the image; the pixels span 0…1, 1…2, …, as in
#   Makie's `image`.
# - `:lowerleftpixel`: the centre of the lower-left pixel; the pixels are centred on 0, 1, …,
#   as in Makie's `heatmap`.
#
# Any other place is given as a fractional index `(j0, i0)`: pixel `(j, i)` has its centre at
# the index `(j, i)` and its corners at `j ± 0.5`, `i ± 0.5`. For example the upper-left corner:

imagedomain((4, 5); origin=(0.5, 0.5))[1, 1]

# ## Units
#
# A domain has no units by default: the coordinates are pixels.

imagedomain((480, 640))

# The pixel size and its unit are given when the domain is made, or declared later with
# [`setunits`](@ref):

dom = imagedomain((480, 640); pixelsize=5.5, units=:um)

#-

units(dom)

# [`convertunits`](@ref) changes the numbers. Between nm, μm, mm, cm and m the factor is an
# exact power of ten, so the coordinates stay round:

convertunits(dom, :mm)

# [`inwavelengths`](@ref) gives coordinates in wavelengths, `x/λ`, with `λ` in the units of
# the domain; the value of `λ` is kept, so the conversion can be undone:

w = inwavelengths(dom, 0.633)

#-

convertunits(w, :um) ≈ dom

# The labels of the axes follow the units:

axislabels(dom), axislabels(w)

# Multiplying a domain by a number drops the units, since a number does not tell what they
# become: dividing by a radius gives dimensionless coordinates, multiplying a grid of spatial
# frequencies by `λ f` gives lengths.

units(dom / 100)

# ## Features in an image
#
# A camera image with two spots, of 20 μm width, and noise. The spots are placed at known points,
# and we find them back.

truth = [(x=-310.0, y=120.0), (x=405.0, y=-230.0)]
spot(p, c) = exp(-((p.x - c.x)^2 + (p.y - c.y)^2) / (2 * 20.0^2))
rng = Xoshiro(1)
img = map(p -> spot(p, truth[1]) + 0.8spot(p, truth[2]), dom) .+ 0.02 .* randn(rng, size(dom))

# The centroid of each spot, in pixel indices, from a window of 31 × 31 pixels around its
# maximum. The indices are fractional:

function centroid(img, J, I; h=15)
    win = CartesianIndices((J-h:J+h, I-h:I+h))
    w = [max(img[k], 0.0) for k in win]
    return (sum(w .* [k[1] for k in win]) / sum(w), sum(w .* [k[2] for k in win]) / sum(w))
end
left = img[:, 1:320]                    # one spot in each half
J1, I1 = Tuple(argmax(left))
J2, I2 = Tuple(argmax(img[:, 321:end])) .+ (0, 320)
c1, c2 = centroid(img, J1, I1), centroid(img, J2, I2)

# [`pointof`](@ref) converts the fractional indices to points of the domain, in μm:

p1, p2 = pointof(dom, c1), pointof(dom, c2)

# The found points agree with the true ones to a fraction of a pixel, and the distance
# between the spots comes out in μm:

hypot((p1 - p2)...), hypot(truth[1].x - truth[2].x, truth[1].y - truth[2].y)

#-

@assert all(abs.(p1 .- (truth[1].x, truth[1].y)) .< 1) && all(abs.(p2 .- (truth[2].x, truth[2].y)) .< 1)

# [`indexof`](@ref) goes the other way, from a point to the fractional index; rounded, it is
# the pixel that contains the point:

round.(Int, indexof(dom, p1)), (J1, I1)

# Circles of radius 50 μm around the found spots, drawn by Makie in the coordinates of the
# domain:

fig = Figure(; size=(700, 520))
ax = Axis(fig[1, 1]; aspect=DataAspect(), axislabels(dom)...)
heatmap!(ax, dom, img; colormap=:grays)
for p in (p1, p2)
    arc!(ax, Point2f(p...), 50, 0, 2π; color=:red)
end
fig
