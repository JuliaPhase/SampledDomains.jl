using SampledDomains
using SampledDomains: SampledDomain, getranges, IMAGEJ_ZERO_YORIGIN, imagej_description
using CoordinateTransformations: AffineMap, LinearMap, Translation
using StaticArrays: SVector, @SMatrix
using AbstractFFTs: fftfreq, fftshift
using TiffImages
using Aqua
using Test

const DATA = joinpath(@__DIR__, "data")

@testset "SampledDomains.jl" begin
    a = CartesianDomain2D(1:3, -5:0.5:-3.5)
    f(x, y) = x * y
    b = SampledDomain(f, a)
    @test b.vals == [
        -5.0 -10.0 -15.0
        -4.5 -9.0 -13.5
        -4.0 -8.0 -12.0
        -3.5 -7.0 -10.5
    ]
    @test b.dom === a
end

@testset "SampledDomain calls f(x, y)" begin
    d = CartesianDomain2D(1:3, 10:10:20)   # 3 samples in x, 2 in y
    fx = SampledDomain((x, y) -> x, d)
    fy = SampledDomain((x, y) -> y, d)
    @test fx.vals == [1 2 3; 1 2 3]
    @test fy.vals == [10 10 10; 20 20 20]
    @test SampledDomain((x, y) -> 100x + y, d).vals ==
        [100x + y for y in d.yrange, x in d.xrange]
    @test_throws DimensionMismatch SampledDomain(zeros(3, 3), d)
    s = SampledDomain(zeros(Float32, 2, 3), d)
    @test s isa SampledDomain{Float32,2,Matrix{Float32},typeof(d)}
end

@testset "CartesianDomain2D: an array of points (x, y)" begin
    d = CartesianDomain2D(1:3, 10:10:20)   # 3 samples in x, 2 in y
    @test d isa AbstractMatrix{SVector{2,Int}}
    @test size(d) == (2, 3)
    @test length(d) == 6
    @test d[2, 3] == SVector(3, 20)        # (x, y)
    @test d[2, 3].x == 3 && d[2, 3].y == 20
    # linear indexing is column-major and agrees with Cartesian indexing and iteration
    for k in 1:length(d)
        @test d[k] == d[Tuple(CartesianIndices(size(d))[k])...]
    end
    @test collect(d) == [SVector(x, y) for y in d.yrange, x in d.xrange]
    @test [p for p in d][:] == [d[k] for k in 1:length(d)]
    @test map(p -> 100p.x + p.y, d) == [100x + y for y in d.yrange, x in d.xrange]
    @test_throws BoundsError d[0]
    @test_throws BoundsError d[7]
    @test_throws BoundsError d[3, 1]

    # concrete types, no allocation
    @test isconcretetype(fieldtype(typeof(d), :xrange))
    @test isconcretetype(fieldtype(typeof(d), :yrange))
    @test isconcretetype(fieldtype(typeof(d), :units))
    dm = make_centered_domain2D(4, 3, 0.5; units=:um)
    @test eltype(dm) == SVector{2,Float64}
    @inferred dm[2, 3]
    getpoint(dom) = @allocated dom[2, 3]
    getpoint(dm)
    @test getpoint(dm) == 0
    @test eltype(CartesianDomain2D(1:3, 0.0:0.5:1.0)) == SVector{2,Float64}

    @test getranges(d) == (1:3, 10:10:20)
    @test d == CartesianDomain2D(1:3, 10:10:20)
    @test d != CartesianDomain2D(1:3, 10:10:30)
    @test d != CartesianDomain2D(1:3, 10:10:20; units=:mm)
end

@testset "show" begin
    d = CartesianDomain2D(1:3, 10:10:20)
    @test repr(d) == "CartesianDomain2D(1:3, 10:10:20)"
    txt = repr(MIME"text/plain"(), d)
    @test startswith(txt, "2×3 CartesianDomain2D{$Int} with points (x, y):")
    @test occursin("x: 1:3", txt) && occursin("y: 10:10:20", txt)
    du = CartesianDomain2D(1:3, 10:10:20; units=:um)
    @test repr(du) == "CartesianDomain2D(1:3, 10:10:20; units=μm)"
    @test occursin("units=μm", repr(MIME"text/plain"(), du))
    dxy = CartesianDomain2D(1:3, 10:10:20; units=(x=:um, y=:mm))
    @test occursin("units=(x=μm, y=mm)", repr(dxy))
    a = AffineDomain2D(d, Translation(1.0, 2.0))
    @test startswith(repr(MIME"text/plain"(), a), "2×3 AffineDomain2D{Float64}")
    @test repr(wavelengthunit(0.633, :um)) == "λ (λ = 633.0 nm)"
end

@testset "imagedomain: orientation of an image, origins" begin
    d = imagedomain((4, 5))                          # 4 rows, 5 columns
    @test size(d) == (4, 5)
    @test step(d.xrange) > 0 && step(d.yrange) < 0   # x right, y up
    @test d[1, 1].x < 0 && d[1, 1].y > 0             # upper left
    @test d[4, 5].x > 0 && d[4, 5].y < 0             # lower right
    @test imagedomain(zeros(4, 5)) == d

    @test collect(d.xrange) == [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test collect(d.yrange) == [1.5, 0.5, -0.5, -1.5]
    f = imagedomain((4, 5); origin=:fftcenter)
    @test f[3, 3] == SVector(0.0, 0.0)               # pixel (rows ÷ 2 + 1, columns ÷ 2 + 1)
    @test collect(f.yrange) == [2.0, 1.0, 0.0, -1.0]
    ll = imagedomain((4, 5); origin=:lowerleft)
    @test ll[4, 1] == SVector(0.5, 0.5)              # the corner of the image is (0, 0)
    llp = imagedomain((4, 5); origin=:lowerleftpixel)
    @test llp[4, 1] == SVector(0.0, 0.0)
    @test llp[1, 5] == SVector(4.0, 3.0)
    t = imagedomain((4, 5); origin=(0.5, 0.5))      # upper-left corner
    @test t[1, 1] == SVector(0.5, -0.5)
    @test_throws ArgumentError imagedomain((4, 5); origin=:middle)
    @test_throws ArgumentError imagedomain((4, 5); pixelsize=-1)

    # pixel size and units
    p = imagedomain((4, 5); pixelsize=(2.0, 3.0), units=:um)
    @test pixelsize(p) == (2.0, 3.0)
    @test collect(p.yrange) == [4.5, 1.5, -1.5, -4.5]
    @test units(p) == (x=LengthUnit(:um), y=LengthUnit(:um))
end

@testset "make_centered_domain2D (y up since 0.2)" begin
    d = make_centered_domain2D(4, 3, 0.5)            # x first
    @test size(d) == (3, 4)
    @test collect(d.xrange) == [-0.75, -0.25, 0.25, 0.75]
    @test collect(d.yrange) == [0.5, 0.0, -0.5]      # descending: row 1 is the top
    @test collect(make_centered_domain2D(4, 3, 0.5, 2.0).yrange) == [2.0, 0.0, -2.0]
    @test d == imagedomain((3, 4); pixelsize=0.5)
    @test units(make_centered_domain2D(4, 3, 0.5; units=:mm)).x == LengthUnit(:mm)
    # the values of a function are those of an image read from a file: row 1 at the top
    up(p) = p.y > 0
    a = map(up, make_centered_domain2D(4, 4, 1.0))
    @test all(a[1:2, :]) && !any(a[3:4, :])
end

@testset "scaling drops the units" begin
    d = make_centered_domain2D(4, 3, 0.5; units=:um)
    s = d * 2
    @test collect(s.xrange) == [-1.5, -0.5, 0.5, 1.5]
    @test collect(s.yrange) == [1.0, 0.0, -1.0]
    @test units(s) == (x=nothing, y=nothing)
    @test (2 * d) == s
    @test (d / 0.5) == s
    @test units(setunits(s, :mm)).y == LengthUnit(:mm)
end

@testset "units: conversion, wavelengths, labels" begin
    @test LengthUnit(:μm) == LengthUnit(:um) == LengthUnit(:micron)
    @test_throws ArgumentError LengthUnit(:inch)
    d = imagedomain((6, 8); pixelsize=(5.5, 4.0), units=:um)
    mm = convertunits(d, :mm)
    @test mm.xrange[1] === -0.01925                 # exact powers of ten
    @test mm.yrange[1] === 0.01
    @test units(mm).x == LengthUnit(:mm)
    @test convertunits(mm, :um) ≈ d
    nm = convertunits(d, :nm)
    @test nm.xrange[1] === -19250.0

    λ = 0.633
    w = inwavelengths(d, λ)
    @test w.xrange ≈ d.xrange ./ λ
    @test units(w).x.name === :λ
    @test units(w).x.metres ≈ 0.633e-6
    @test convertunits(w, :um).xrange ≈ d.xrange
    @test axislabels(w) == (xlabel="x, λ", ylabel="y, λ")
    @test axislabels(d) == (xlabel="x, μm", ylabel="y, μm")
    @test axislabels(imagedomain((2, 2))) == (xlabel="x", ylabel="y")
    @test axislabels(setunits(d, (x=:mm, y=nothing))) == (xlabel="x, mm", ylabel="y")

    @test_throws ArgumentError convertunits(imagedomain((2, 2)), :mm)   # no units
    @test_throws ArgumentError inwavelengths(imagedomain((2, 2)), λ)
    @test_throws ArgumentError convertunits(dualDomain(d), :mm)          # μm⁻¹ to mm
    @test convertunits(dualDomain(d), inv(LengthUnit(:mm))).xrange ≈ dualDomain(d).xrange .* 1000
end

@testset "pointof, indexof" begin
    d = imagedomain((6, 8); pixelsize=(5.5, 4.0), origin=:lowerleft, units=:um)
    for I in CartesianIndices(d)
        @test pointof(d, I) ≈ d[I]
        @test all(indexof(d, d[I]) .≈ Tuple(I))
    end
    @test pointof(d, 6.5, 0.5) ≈ SVector(0.0, 0.0)   # the corner of the image
    @test pointof(d, (1.25, 2.75)) ≈ SVector(2.25 * 5.5, 5.25 * 4.0)
    @test round.(Int, indexof(d, SVector(12.0, 7.0))) == (5, 3)

    grid = imagedomain((5, 7); pixelsize=0.5)
    rot = LinearMap(@SMatrix [0.0 -1.0; 1.0 0.0])
    a = AffineDomain2D(grid, Translation(0.3, -0.1) ∘ rot)
    @test all(pointof(a, Tuple(I)...) ≈ a[I] for I in CartesianIndices(a))
    @test all(all(indexof(a, a[I]) .≈ Tuple(I)) for I in CartesianIndices(a))
end

@testset "AffineDomain2D" begin
    grid = make_centered_domain2D(4, 3, 0.5)
    rot90 = LinearMap(@SMatrix [0.0 -1.0; 1.0 0.0])
    t = Translation(0.1, -0.2)

    dr = AffineDomain2D(grid, rot90)
    @test size(dr) == size(grid)
    @test eltype(dr) == SVector{2,Float64}
    for I in CartesianIndices(grid)
        p = grid[I]
        @test dr[I] == SVector(-p.y, p.x)
    end

    # composition: rotate, then shift; the map is stored as an AffineMap
    da = AffineDomain2D(grid, t ∘ rot90)
    @test da.map isa AffineMap
    @test all(da[I] ≈ rot90(grid[I]) + SVector(0.1, -0.2) for I in CartesianIndices(grid))
    # the inverse map gives back the grid points
    @test all(inv(da.map)(da[I]) ≈ grid[I] for I in CartesianIndices(grid))
    # a translation alone; a plain matrix is converted to a static one
    @test AffineDomain2D(grid, t)[1, 1] ≈ grid[1, 1] + SVector(0.1, -0.2)
    dm = AffineDomain2D(grid, AffineMap([1.0 0.0; 0.0 -1.0], [0.0, 0.0]))   # reflection in y
    @test dm[1, 2] == SVector(grid[1, 2].x, -grid[1, 2].y)
    @test isbitstype(typeof(dm.map))
    getpoint(dom) = @allocated dom[2, 3]
    getpoint(dm)
    @test getpoint(dm) == 0
    @inferred dm[2, 3]

    # scaling acts on the physical coordinates
    @test all((da * 2)[I] ≈ 2da[I] for I in CartesianIndices(grid))
    @test 2 * da == da * 2

    # functions on an affine domain see the physical coordinates
    s = SampledDomain((x, y) -> x, dr)
    @test s.vals == [-p.y for p in grid]

    # units: of the physical frame, one for both axes
    gu = setunits(grid, :mm)
    au = AffineDomain2D(gu, rot90)
    @test units(au).x == LengthUnit(:mm)
    @test units(au * 2) == (x=nothing, y=nothing)
    @test convertunits(au, :um)[2, 3] ≈ 1000 * au[2, 3]
    @test_throws ArgumentError AffineDomain2D(setunits(grid, (x=:mm, y=:um)), rot90)
    @test AffineDomain2D(grid, rot90; units=:um).units.x == LengthUnit(:um)
end

@testset "dualRange, dualDomain" begin
    # even and odd lengths, upsampling; the order is that of fftshift
    @test collect(dualRange(0:0.5:1.5)) == [-1.0, -0.5, 0.0, 0.5]
    @test collect(dualRange(0:0.5:2)) == [-0.8, -0.4, 0.0, 0.4, 0.8]
    @test collect(dualRange(0:0.5:1.5, 2)) == [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    dd = dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), (1, 2))
    @test collect(dd.xrange) == [-1.0, -0.5, 0.0, 0.5]
    @test collect(dd.yrange) ≈ [-4, -8 / 3, -4 / 3, 0, 4 / 3, 8 / 3]
    # a descending range (y up) gives a descending dual range, zero at n ÷ 2 + 1
    r = dualRange(1.5:-0.5:0)
    @test r isa AbstractRange
    @test collect(r) == [1.0, 0.5, 0.0, -0.5]
    @test units(dualDomain(make_centered_domain2D(4, 4, 1.0; units=:um))).x ==
        inv(LengthUnit(:um))
end

@testset "dualDomain: orientation of the FFT of an array on an image domain" begin
    # a plane wave exp(2πi (u x + v y)) with v > 0 (phase growing upward) has its spectrum at
    # (u, v) of the dual domain, i.e. in the upper half for a domain with y up. Plain DFT.
    n = 8
    d = make_centered_domain2D(n, n, 0.5)
    dd = dualDomain(d)
    u, v = dd.xrange[n ÷ 2 + 2], dd.yrange[2]        # a positive u and the largest v
    @test u > 0 && v > 0
    a = map(p -> cispi(2 * (u * p.x + v * p.y)), d)
    dft = [sum(a[j, i] * cispi(-2 * (q.x * d[j, i].x + q.y * d[j, i].y)) for j in 1:n, i in 1:n)
           for q in dd]
    jmax, imax = Tuple(argmax(abs.(dft)))
    @test dd[jmax, imax] ≈ SVector(u, v)
    @test jmax == 2                                   # upper rows: y of the spectrum up
end

@testset "dualRange: fftfreq, accurate to one ulp" begin
    for (n, s, q) in ((8, 0.1, 1), (9, 0.1, 1), (8, 0.013, 2), (7, 2.5e-3, 3), (101, 0.37, 2))
        x = range(0; step=s, length=n)
        r = dualRange(x, q)
        @test r isa AbstractRange
        @test r == fftshift(fftfreq(n * q, q / s))
        ks = (-((n * q) ÷ 2)):(cld(n * q, 2) - 1)
        exact = [BigFloat(k) / (BigFloat(s) * n) for k in ks]
        @test all(abs(r[m] - exact[m]) <= eps(r[m]) for m in eachindex(ks) if !iszero(ks[m]))
        @test r[(n * q) ÷ 2 + 1] === 0.0
    end
    @test dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), 2) ==
        dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), (2, 2))
end

@testset "ImageJ calibration" begin
    # the four origins of an image domain, y up
    for (o, xo, yo) in ((:center, 4.0, 3.0), (:fftcenter, 4.5, 3.5), (:lowerleft, 0.0, 6.0),
        (:lowerleftpixel, 0.5, 5.5))
        cal = imagej_calibration(imagedomain((6, 8); pixelsize=(5.5, 4.0), origin=o, units=:um))
        @test cal.xorigin ≈ xo atol = 1e-12
        @test cal.yorigin ≈ yo atol = 1e-12
        @test cal.inverty && cal.unit == "micron" && cal.yunit === nothing
        @test cal.pixelwidth == 5.5 && cal.pixelheight == 4.0
    end
    # the quirk of ImageJ: an origin on the upper edge with y up
    q = imagedomain((6, 8); origin=(0.5, 4.5), units=:um)
    @test imagej_calibration(q).yorigin == IMAGEJ_ZERO_YORIGIN
    @test imagej_domain((6, 8), imagej_description(imagej_calibration(q)), 1.0, 1.0) ≈ q
    # y down (ascending yrange) is written without inverty
    down = CartesianDomain2D(1.0:1.0:8.0, 1.0:1.0:6.0; units=:mm)
    cd = imagej_calibration(down)
    @test !cd.inverty && cd.unit == "mm"
    @test imagej_domain((6, 8), imagej_description(cd), 1.0, 1.0) ≈ down
    # different units, wavelengths
    @test imagej_calibration(setunits(down, (x=:mm, y=:um))).yunit == "micron"
    w = inwavelengths(imagedomain((6, 8); units=:um), 0.5)
    cw = imagej_calibration(w)
    @test cw.unit == "lambda" && cw.wavelength ≈ 0.5e-6
    @test units(imagej_domain((6, 8), imagej_description(cw), cw.pixelwidth, cw.pixelheight)) ==
        units(w)
    @test_throws ArgumentError imagej_calibration(CartesianDomain2D(2:-1:1, 1:2))
end

@testset "ImageJ files" begin
    # files written by ImageJ 1.54g (test/data/README.md): value = x + 10y of the 0-based
    # pixel (x, y) of ImageJ, pixel 5.5 × 4 micron, origin (4, 3) pixels. ImageJ measured the
    # centres of the pixels (3, 1) (value 13) and (6, 4) (value 46).
    for (file, c13, c46, ydir) in (("cal_inverty.tif", (-2.75, 6.0), (13.75, -6.0), -1),
        ("cal_plain.tif", (-2.75, -6.0), (13.75, 6.0), 1))
        img, d = load_imagej(joinpath(DATA, file))
        @test size(img) == size(d) == (6, 8)
        @test units(d).x == LengthUnit(:um)
        @test sign(step(d.yrange)) == ydir
        @test Float64(img[2, 4]) == 13 && Float64(img[5, 7]) == 46
        @test all(Float64(img[j, i]) == (i - 1) + 10(j - 1) for j in 1:6, i in 1:8)
        @test d[2, 4] ≈ SVector(c13) rtol = 1e-5     # ImageJ writes 5.5 as 1/0.181818
        @test d[5, 7] ≈ SVector(c46) rtol = 1e-5
        @test imagej_domain(joinpath(DATA, file)) == d
    end

    # round trip through a file written by us
    mktempdir() do dir
        dom = imagedomain((6, 8); pixelsize=(5.5, 4.0), origin=:lowerleftpixel, units=:um)
        a = reshape(collect(1.0:48.0), 6, 8)
        f = joinpath(dir, "a.tif")
        save_imagej(f, a, dom)
        img, d2 = load_imagej(f)
        @test d2 == dom                                # the pixel size is an exact rational
        @test Float64.(img) ≈ a
        @test imagej_domain(f) == dom
        # an origin on the upper edge
        q = imagedomain((6, 8); origin=(0.5, 4.5), units=:um)
        save_imagej(f, a, q)
        @test imagej_domain(f) ≈ q
        @test_throws DimensionMismatch save_imagej(f, zeros(3, 3), q)
        # a file without the calibration of ImageJ
        g = joinpath(dir, "plain.tif")
        TiffImages.save(g, TiffImages.DenseTaggedImage(TiffImages.Gray{Float32}.(rand(Float32, 4, 5))))
        @test imagej_domain(g) == imagedomain((4, 5))
    end
end

@testset "Aqua" begin
    Aqua.test_all(SampledDomains)
end
