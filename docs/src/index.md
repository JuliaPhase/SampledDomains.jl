```@meta
CurrentModule = SampledDomains
```

# SampledDomains

SampledDomains attaches physical coordinates to array elements.

## Array axis convention

Arrays defined on a [`CartesianDomain2D`](@ref) follow the Julia image convention:
the **first** array index runs along `y` and the **second** along `x`.

```
a[j, i] = f(x[i], y[j])        size(a) == (length(yrange), length(xrange))
```

Row 1 corresponds to `yrange[1]`. This is the same index order as for images loaded
with `FileIO`/`Images` and the layout expected by FFT-based code, so arrays on a domain can be
combined with image data without transposition.

```julia
using SampledDomains

dom = CartesianDomain2D(-1:0.1:1, -0.5:0.1:0.5)   # 21 samples in x, 11 in y
a = [x + 2y for y in dom.yrange, x in dom.xrange]  # a[j, i] = x[i] + 2y[j]
size(a)                                            # (11, 21)
```

### Displaying arrays

Makie's `heatmap(A)` puts the **first** index on the horizontal axis, so it shows
`a` with `x` and `y` interchanged. To get `x` to the right and `y` up, transpose:

```julia
heatmap(dom.xrange, dom.yrange, a')    # axis values are shown, y increases upward
```

With `PhasePlots` this is done for you: `showarray(dom, a)` shows `a` with `x` to the right and `y` up,
real coordinates on the axes, and axes labelled `x` and `y`.

### Evaluating a function on a domain

`SampledDomain(f, dom)` calls `f(x, y)` and stores the result with the layout above:

```julia
s = SampledDomain((x, y) -> x + 2y, dom)
s.vals[j, i] == dom.xrange[i] + 2dom.yrange[j]
```

### Relation to images read from files

An image loaded from a file has row 1 at the **top**, so a Cartesian frame with `y` up
is `y_user = -ind1`. An array on a domain with an ascending `yrange` has row 1 at the
**bottom** of a y-up plot. The two therefore differ by a vertical flip. To work in the
image frame, either use a descending `yrange`, or flip the array (`reverse(a; dims=1)`)
before comparing. Both frames use the same index order, so no transposition is needed.

## API

```@index
```

```@autodocs
Modules = [SampledDomains]
```
