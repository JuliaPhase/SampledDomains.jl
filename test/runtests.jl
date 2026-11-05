using SampledDomains
using SampledDomains: SampledDomain, getranges
using CoordinateTransformations: AffineMap, LinearMap, Translation
using StaticArrays: SVector, @SMatrix
using AbstractFFTs: fftfreq, fftshift
using Aqua
using Test

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
    dm = make_centered_domain2D(4, 3, 0.5)
    @test eltype(dm) == SVector{2,Float64}
    @inferred dm[2, 3]
    getpoint(dom) = @allocated dom[2, 3]
    getpoint(dm)
    @test getpoint(dm) == 0
    @test eltype(CartesianDomain2D(1:3, 0.0:0.5:1.0)) == SVector{2,Float64}

    @test getranges(d) == (1:3, 10:10:20)
    @test d == CartesianDomain2D(1:3, 10:10:20)
    @test d != CartesianDomain2D(1:3, 10:10:30)
end

@testset "show" begin
    d = CartesianDomain2D(1:3, 10:10:20)
    @test repr(d) == "CartesianDomain2D(1:3, 10:10:20)"
    txt = repr(MIME"text/plain"(), d)
    @test startswith(txt, "2×3 CartesianDomain2D{$Int} with points (x, y):")
    @test occursin("x: 1:3", txt) && occursin("y: 10:10:20", txt)
    a = AffineDomain2D(d, Translation(1.0, 2.0))
    @test startswith(repr(MIME"text/plain"(), a), "2×3 AffineDomain2D{Float64}")
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
end

@testset "centred domains, scaling" begin
    d = make_centered_domain2D(4, 3, 0.5)
    @test collect(d.xrange) == [-0.75, -0.25, 0.25, 0.75]
    @test collect(d.yrange) == [-0.5, 0.0, 0.5]
    @test collect(make_centered_domain2D(4, 3, 0.5, 2.0).yrange) == [-2.0, 0.0, 2.0]
    s = d * 2
    @test collect(s.xrange) == [-1.5, -0.5, 0.5, 1.5]
    @test collect(s.yrange) == [-1.0, 0.0, 1.0]
    @test (2 * d).xrange == s.xrange && (2 * d).yrange == s.yrange
end

@testset "dualRange, dualDomain" begin
    # even and odd lengths, upsampling; the order is that of fftshift
    @test collect(dualRange(0:0.5:1.5)) == [-1.0, -0.5, 0.0, 0.5]
    @test collect(dualRange(0:0.5:2)) == [-0.8, -0.4, 0.0, 0.4, 0.8]
    @test collect(dualRange(0:0.5:1.5, 2)) == [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    dd = dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), (1, 2))
    @test collect(dd.xrange) == [-1.0, -0.5, 0.0, 0.5]
    @test collect(dd.yrange) ≈ [-4, -8 / 3, -4 / 3, 0, 4 / 3, 8 / 3]
end

@testset "dualRange agrees with fftfreq" begin
    for (n, s, q) in ((8, 0.1, 1), (9, 0.1, 1), (8, 0.013, 2), (7, 2.5e-3, 3), (101, 0.37, 2))
        x = range(0; step=s, length=n)
        @test dualRange(x, q) ≈ fftshift(fftfreq(n * q, q / s)) rtol = 1e-12
    end
    @test dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), 2) ==
        dualDomain(CartesianDomain2D(0:0.5:1.5, 0:0.25:0.5), (2, 2))
end

@testset "Aqua" begin
    Aqua.test_all(SampledDomains)
end
