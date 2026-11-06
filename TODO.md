# SampledDomains — TODO

`AffineDomain2D` is implemented in 0.2.0 (on CoordinateTransformations.jl). Still open:

## `dualDomain` of an affine domain

The dual of the linear part `A` is `(A⁻¹)ᵀ`; the translation does not shift the dual grid (it
becomes a phase factor). To be implemented when it is needed (rotated cameras in
PhaseRetrieval), with a test against the FFT of a rotated grid.

## Makie: axis labels from the domain

The extension draws `heatmap(dom, a)` etc. with the axes of the domain; the labels with units
are passed by hand, `Axis(...; axislabels(dom)...)`. Setting them automatically would need a
plot recipe that owns the axis. `AffineDomain2D` is not drawn yet (`surface` with the matrices
of `x` and `y` could do it).

## DimensionalData extension (if wanted)

Conversion of `(a, dom)` to a `DimArray` with `X` and `Y` dimensions, for plots with physical
axes or for exchange with other software.

## Users of the new interface

After the change of orientation (0.2): check every user against results computed with `v0.1.0`
(design: `Phase/design/SampledDomains.md` §8.1).


- PhaseBases: `coordmap` → `AffineDomain2D` (the `CoordMap` item of its `TODO.md`).
- PhaseRetrieval: `get_polarization_magnitudes` (VectorialPSF.jl) could iterate over the points
  of the domain.
- DeformableMirrors: the beam on the mirror as an `AffineDomain2D`.
