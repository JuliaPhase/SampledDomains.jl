module SampledDomainsMakieExt

# Plots of an array on a domain, with the axes of the domain. Makie expects a[ix, iy] (the
# first index horizontal), while arrays on a domain are a[j, i] (row j = y, column i = x), so
# the array is transposed here; a descending yrange (y up) is drawn correctly by Makie.
# Without this, heatmap(x, y, a) draws a wrong picture without an error.

using SampledDomains
using Makie

const _GridPlots = Union{Type{<:Heatmap},Type{<:Contour},Type{<:Contourf},Type{<:Surface}}

Makie.convert_arguments(P::_GridPlots, dom::CartesianDomain2D, a::AbstractMatrix) =
    Makie.convert_arguments(P, dom.xrange, dom.yrange, permutedims(a))

function Makie.convert_arguments(P::Type{<:Image}, dom::CartesianDomain2D, a::AbstractMatrix)
    # image takes the outer edges; the half pixels around the first and last centres
    x, y = dom.xrange, dom.yrange
    hx, hy = step(x) / 2, step(y) / 2
    return Makie.convert_arguments(P, (first(x) - hx, last(x) + hx),
        (first(y) - hy, last(y) + hy), permutedims(a))
end

end
