using SampledDomains
using Test

@testset "SampledDomains.jl" begin
    a= SampledDomains.CartesianDomain2D(1:3, -5:.5:-3.5)
    f(x,y) = x*y
    b = SampledDomains.SampledDomain(f, a)
    @test b.vals ==[
        -5.0  -10.0  -15.0
        -4.5   -9.0  -13.5
        -4.0   -8.0  -12.0
        -3.5   -7.0  -10.5
    ] 
    @test b.dom === a
end

@testset "SampledDomain calls f(x, y)" begin
    d = SampledDomains.CartesianDomain2D(1:3, 10:10:20)   # 3 samples in x, 2 in y
    fx = SampledDomains.SampledDomain((x, y) -> x, d)
    fy = SampledDomains.SampledDomain((x, y) -> y, d)
    @test fx.vals == [1 2 3; 1 2 3]
    @test fy.vals == [10 10 10; 20 20 20]
    @test SampledDomains.SampledDomain((x, y) -> 100x + y, d).vals ==
        [100x + y for y in d.yrange, x in d.xrange]
    @test_throws ErrorException SampledDomains.SampledDomain(zeros(3, 3), d)
end

@testset "CartesianDomain2D indexing" begin
    d = SampledDomains.CartesianDomain2D(1:3, 10:10:20)   # 3 samples in x, 2 in y
    @test size(d) == (2, 3)
    @test length(d) == 6
    @test d[2, 3] == [20, 3]          # [y, x]
    # linear indexing is column-major and agrees with Cartesian indexing and iteration
    for k in 1:length(d)
        @test d[k] == d[Tuple(CartesianIndices(size(d))[k])...]
    end
    @test [d[k] for k in 1:length(d)] == [[y, x] for (y, x) in d][:]
    @test_throws BoundsError d[0]
    @test_throws BoundsError d[7]
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
