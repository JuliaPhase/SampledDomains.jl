# ```@meta
# CurrentModule = SampledDomains
# ```

using SampledDomains
using CairoMakie
CairoMakie.activate!(; type="png")

# # Arrays and pictures
#
# A two-dimensional array can be shown in three ways, and they disagree. This page shows the
# same array in all three, so that the convention of SampledDomains can be seen rather than
# remembered.
#
# ## A test object
#
# The test object is the letter F: it looks different after every flip, rotation or
# transposition. We sample it on a grid of 120 points in `x` and 80 points in `y`, with step
# 0.05, centred on zero.

dom = make_centered_domain2D(120, 80, 0.05)

# A domain is the array of the points of the samples. `dom[j, i]` is the point `(x, y)` of the
# sample `[j, i]`; the first index runs along `y`, the second along `x`, so the domain has 80
# rows and 120 columns:

size(dom)

#-

dom[1, 1]

# The first sample lies at the lowest `x` and the lowest `y`. The coordinates can be read by
# name:

p = dom[1, 120]
p.x, p.y

# The letter is described in `x` and `y`, and `map` evaluates it at every point of the domain.
# The result is an ordinary array of the size of the domain.

letterF(p) =
    (-1.5 <= p.x <= -1.0 && -1.5 <= p.y <= 1.5) ||    # stem
    (-1.5 <= p.x <= 1.5 && 1.0 <= p.y <= 1.5) ||      # top bar
    (-1.5 <= p.x <= 0.5 && 0.0 <= p.y <= 0.5)         # middle bar

a = Float64.(map(letterF, dom))
size(a)

# ## Three pictures of one array

x, y = dom.xrange, dom.yrange
fig = Figure(; size=(1000, 330))

ax1 = Axis(fig[1, 1]; title="heatmap(a)", xlabel="first index j", ylabel="second index i",
    aspect=DataAspect())
heatmap!(ax1, a; colormap=:grays)

ax2 = Axis(fig[1, 2]; title="heatmap(x, y, a')", xlabel="x", ylabel="y", aspect=DataAspect())
heatmap!(ax2, x, y, a'; colormap=:grays)
scatter!(ax2, [dom[1, 1].x], [dom[1, 1].y]; color=:red, markersize=12)
text!(ax2, dom[1, 1].x, dom[1, 1].y; text=" a[1, 1]", color=:red, align=(:left, :bottom))

ax3 = Axis(fig[1, 3]; title="image viewer: row 1 on top", xlabel="x", ylabel="y",
    aspect=DataAspect(), yreversed=true)
heatmap!(ax3, x, y, a'; colormap=:grays)
scatter!(ax3, [dom[1, 1].x], [dom[1, 1].y]; color=:red, markersize=12)

fig

# - Left: Makie puts the first index of the array on the horizontal axis. The picture is the
#   letter transposed: `x` and `y` are interchanged.
# - Middle: transposing the array and giving the coordinates shows the letter as it was
#   defined, with `x` to the right and `y` up. The red point is `a[1, 1]`, at the lowest `x` and
#   `y`.
# - Right: a file viewer (and JuliaImages) draws row 1 at the top. The letter is upside down:
#   the difference from the middle picture is a vertical flip, not a transposition.
#
# The middle picture is the one to use for physics. The PhasePlots function
# `showarray(dom, a)` draws it directly.
#
# The right picture explains a frequent sign error. An array read from an image file has row 1
# at the top, so its `y` axis points down. To bring it into the frame of the middle picture,
# flip it, `reverse(img; dims=1)`, or describe it by a domain with a descending `yrange`.

# ## The array in the REPL
#
# The REPL is a fourth way to look at an array. On a coarse grid of 17 × 17 points with step 0.25
# the letter is small enough to be printed:

small = CartesianDomain2D(-2:0.25:2, -2:0.25:2)
F = Int.(map(letterF, small))

# The REPL prints row 1 on top, like the image viewer of the right picture above, so the
# letter is upside down. Row 1 is `y = -2`, the lowest `y`. Flipping the rows gives the letter
# as it was defined, with `y` up:

reverse(F; dims=1)

# ## Stacks of arrays
#
# For more than two dimensions the spatial axes come first, in the order `(y, x)`, and the
# others follow, as in JuliaImages. The domain describes the first two axes, and every slice is
# an array on it. Here are four frames of a growing disc:

radii = [0.5, 1.0, 1.5, 2.0]
stack = zeros(size(dom)..., length(radii))
for (frame, r) in zip(eachslice(stack; dims=3), radii)
    frame .= map(p -> hypot(p.x - 0.5, p.y) <= r, dom)
end
size(stack)

#-

fig = Figure(; size=(1000, 220))
for (k, frame) in enumerate(eachslice(stack; dims=3))
    ax = Axis(fig[1, k]; title="frame $k, r = $(radii[k])", aspect=DataAspect())
    heatmap!(ax, x, y, frame'; colormap=:grays)
    hidedecorations!(ax)
end
fig
