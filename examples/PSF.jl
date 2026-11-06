# ```@meta
# CurrentModule = SampledDomains
# ```

using SampledDomains
using CairoMakie
using FFTW
CairoMakie.activate!(; type="png")

# # PSF in wavelengths
#
# The PSF of a lens is the squared modulus of the Fourier transform of the field in its pupil.
# This page computes it with an FFT, finds the coordinates of the focal plane from the
# Fourier-dual domain, and gives them in wavelengths, `x/λ`: the numbers become dimensionless
# and of the order of the PSF size. Two checks: the first dark ring of the Airy pattern, and
# the direction in which a tilt of the wavefront moves the PSF.
#
# ## The pupil
#
# A lens of diameter `D = 10` mm and focal length `f = 200` mm, light of `λ = 0.633` μm. The
# pupil is sampled with 64 pixels across the diameter, on a grid 8 times larger than the pupil
# (the zero padding sets the sampling of the focal plane). The origin is the pixel of the zero
# frequency of `fftshift`, so that the centre of the pupil goes to the first element of the
# array with `ifftshift`.

D, f, λ = 10.0, 200.0, 0.633e-3          # mm
n = 512
pupil = imagedomain((n, n); pixelsize=D / 64, origin=:fftcenter, units=:mm)

#-

aperture = map(p -> hypot(p.x, p.y) <= D / 2, pupil)
psf(field) = abs2.(fftshift(fft(ifftshift(field))))
I0 = psf(complex.(aperture))
I0 ./= maximum(I0);

# ## The focal plane
#
# The `fftshift`ed FFT of an array on `pupil` is an array on the dual domain, the grid of the
# spatial frequencies, in mm⁻¹. Its orientation is that of the pupil: `y` up.

ν = dualDomain(pupil)

# A spatial frequency `ν` of the pupil field arrives in the focal plane at `x = λ f ν`.
# Multiplying by a number drops the units, so they are declared again; then the coordinates
# are converted to wavelengths:

focal = setunits(ν * (λ * f), :mm)
focalλ = inwavelengths(focal, λ)

# In wavelengths, the first dark ring of the Airy pattern has the radius `1.22 f/D`:

r1 = 1.22 * f / D

# ## A tilt of the wavefront
#
# A wavefront tilted by a small angle `α` upward, `W(x, y) = α y` (the optical path grows with
# `y`), is a plane wave going up; the field in the pupil is `exp(2πi W/λ)`. The PSF moves to
# `y = f α` in the focal plane. With `α = 1.5 × 10⁻⁴` this is 30 μm:

α = 1.5e-4
W = map(p -> α * p.y, pupil)
I1 = psf(aperture .* cis.(2π .* W ./ λ))
I1 ./= maximum(I1)
expected = f * α / λ                      # in wavelengths

#-

fig = Figure(; size=(1000, 420))
lims = (-60, 60, -60, 60)
log10I(I) = log10.(I .+ 1e-5)
ax1 = Axis(fig[1, 1]; title="PSF, log scale", aspect=DataAspect(), limits=lims,
    axislabels(focalλ)...)
heatmap!(ax1, focalλ, log10I(I0); colormap=:inferno, colorrange=(-5, 0))
arc!(ax1, Point2f(0, 0), r1, 0, 2π; color=:cyan, linestyle=:dash)
text!(ax1, r1 * 0.75, r1 * 0.75; text="1.22 f/D", color=:cyan, align=(:left, :bottom))
ax2 = Axis(fig[1, 2]; title="tilted wavefront, W = α y", aspect=DataAspect(), limits=lims,
    axislabels(focalλ)...)
heatmap!(ax2, focalλ, log10I(I1); colormap=:inferno, colorrange=(-5, 0))
scatter!(ax2, [Point2f(0, expected)]; color=:cyan, marker=:cross, markersize=14)
text!(ax2, 3, expected; text="y = f α", color=:cyan, align=(:left, :center))
fig

# The dashed circle at `1.22 f/D` lies on the first dark ring, and the tilted PSF is centred
# at `y = f α`, above the axis, as the plane wave goes up. Both are checked numerically: the
# minimum of the intensity along the positive `x` axis near the ring, and the centroid of the
# tilted PSF.

row0 = findfirst(==(0.0), focalλ.yrange)                  # the row y = 0
xs = focalλ.xrange
near = findall(x -> 0.8r1 < x < 1.2r1, xs)
xmin = xs[near[argmin(I0[row0, near])]]
centroid_y = sum(p.y * I for (p, I) in zip(focalλ, I1)) / sum(I1)
(xmin, r1), (centroid_y, expected)

#-

@assert abs(xmin - r1) < 2step(xs) && abs(centroid_y - expected) < 0.5

# The direction of the shift depends on the sign conventions of optics (here the field
# `exp(2πi W/λ)`, a wave `exp(i(kz - ωt))`, and the forward FFT for the propagation to the
# focal plane); the domains only make sure that `y` up in the pupil and `y` up in the focal
# plane are the same direction as in the images of the setup.
