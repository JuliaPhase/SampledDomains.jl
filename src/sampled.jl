"""
    SampledDomain(vals::AbstractArray, dom::AbstractDomain)
    SampledDomain(f::Function, dom::AbstractDomain)

Values of a function on a domain: the fields `vals` and `dom`, with `size(vals) == size(dom)`.

When constructed from a function, `f` is called with the coordinates of each point in their
natural order, `f(x, y)` for a two-dimensional domain, and `vals` has the layout of the domain:
`vals[j, i] = f(dom[j, i]...)`, i.e. `f(x[i], y[j])` for a [`CartesianDomain2D`](@ref).
"""
struct SampledDomain{T,N,A<:AbstractArray{T,N},D<:AbstractDomain}
    vals::A
    dom::D
    function SampledDomain(vals::AbstractArray, dom::AbstractDomain)
        size(vals) == size(dom) || throw(
            DimensionMismatch(
                "values of size $(size(vals)) on a domain of size $(size(dom))"
            ),
        )
        return new{eltype(vals),ndims(vals),typeof(vals),typeof(dom)}(vals, dom)
    end
end

SampledDomain(f::Function, dom::AbstractDomain) = SampledDomain(map(p -> f(p...), dom), dom)
