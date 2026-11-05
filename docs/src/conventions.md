```@meta
CurrentModule = SampledDomains
```

# [Conventions](@id conventions)

## Index order

Arrays on a [`CartesianDomain2D`](@ref) follow the image convention: the **first** array index
runs along `y` and the **second** along `x`.

```
a[j, i] = f(x[i], y[j])        size(a) == (length(yrange), length(xrange))
dom[j, i] == SVector(x[i], y[j])
```

The index order is that of images loaded with FileIO/Images and of FFT-based code, so arrays
on a domain combine with image data without transposition. The coordinates of a point are
always in the natural order `(x, y)`, and they can be read by name: `p.x`, `p.y`.

## Direction of `y`

Row 1 corresponds to `yrange[1]`. With an ascending `yrange`, row 1 is the lowest `y`.

An image loaded from a file has row 1 at the **top**, so in a frame with `y` up its rows run
downwards. The two therefore differ by a vertical flip, not by a transposition. To work in the
frame of an image, either use a descending `yrange`, or flip the array (`reverse(a; dims=1)`)
before comparing. This is the usual source of a wrong sign of a vertical tilt or decentre.

## More than two dimensions

The spatial axes come first, in the order `(y, x)`, followed by the other axes (frame,
defocus, mode number), as in JuliaImages. The domain describes the first two axes:

```julia
dom = make_centered_domain2D(64, 64, 0.1)
stack = zeros(size(dom)..., 5)          # five frames on dom
for frame in eachslice(stack; dims=3)   # each slice is an array on dom
    frame .= map(p -> p.x^2 + p.y^2, dom)
end
```

The order of the axes is fixed, not chosen by the user: every function that takes a domain
would otherwise have to handle both orders. Conversions happen at the boundaries, for example
when plotting.

## Plotting

Makie's `heatmap(A)` puts the **first** index on the horizontal axis, so it shows `a` with `x`
and `y` interchanged. To get `x` to the right and `y` up, transpose:

```julia
heatmap(dom.xrange, dom.yrange, a')    # axis values are shown, y increases upward
```

With PhasePlots this is done for you: `showarray(dom, a)` shows `a` with `x` to the right and
`y` up, the coordinates on the axes, and the axes labelled `x` and `y`.
