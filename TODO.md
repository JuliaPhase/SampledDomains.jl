# SampledDomains — TODO

## AffineDomain2D — arbitrary linear coordinate transforms

**Motivation:** `CartesianDomain2D` stores axis-aligned `xrange`/`yrange`. When the
physical coordinate frame is rotated or reflected relative to the array axes (e.g.
TIFF files from a camera with a rotated sensor), there is no way to express this in
the domain. Downstream code (Zernike basis in PhaseBases, polarization magnitudes in
PhaseRetrieval) must work around it with ad-hoc `coordmap` arguments.

**Proposed type:**

```julia
struct AffineDomain2D{TA<:AbstractMatrix, TC}
    origin::TC          # (x0, y0) — physical origin of array element [1,1]
    axes::TA            # 2×2 matrix: columns are physical directions of array axes
    size::NTuple{2,Int} # (nrows, ncols)
end
```

`dom[i, j]` returns `origin + axes * [j-1, i-1]` (or similar convention), giving
physical `(x, y)` for every array index.

**Requirements:**
- `*(scalar, dom)` — scale physical coordinates
- `dualDomain(dom, q)` — Fourier-dual of an affine domain
- Iteration protocol consistent with `CartesianDomain2D`
- `make_centered_affine_domain2D(ims, rotation_angle)` constructor

**Impact if implemented:**
- `ZernikeBW(dom::AffineDomain2D, d, order)` would pick up the rotation automatically;
  no `coordmap` keyword needed.
- `get_polarization_magnitudes` in PhaseRetrieval/VectorialPSF.jl would iterate over
  physically correct `(x, y)` pairs without any changes.
- The `coordmap` workaround in PhaseBases could be deprecated.

**Current workaround:** `coordmap` keyword in `ZernikeBW` (PhaseBases), used via
`SimParams.coordmap` in Feedback14AMI.
