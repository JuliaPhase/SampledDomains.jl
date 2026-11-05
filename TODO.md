# SampledDomains — TODO

`AffineDomain2D` is implemented in 0.2.0 (on CoordinateTransformations.jl). Still open:

## `dualDomain` of an affine domain

The dual of the linear part `A` is `(A⁻¹)ᵀ`; the translation does not shift the dual grid (it
becomes a phase factor). To be implemented when it is needed (rotated cameras in
PhaseRetrieval), with a test against the FFT of a rotated grid.

## Makie recipe (package extension)

`heatmap(dom, a)` that transposes the array, puts the coordinates on the axes and labels them
`x` and `y`, so that the order of indices of the image convention never reaches the user.
PhasePlots `showarray(dom, a)` does this now.

## DimensionalData extension (if wanted)

Conversion of `(a, dom)` to a `DimArray` with `X` and `Y` dimensions, for plots with physical
axes or for exchange with other software.

## Users of the new interface

- PhaseBases: `coordmap` → `AffineDomain2D` (the `CoordMap` item of its `TODO.md`).
- PhaseRetrieval: `get_polarization_magnitudes` (VectorialPSF.jl) could iterate over the points
  of the domain.
- DeformableMirrors: the beam on the mirror as an `AffineDomain2D`.
