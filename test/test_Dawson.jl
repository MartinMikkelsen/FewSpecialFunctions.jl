using FewSpecialFunctions
import FewSpecialFunctions: dawson
using Test

@testset "Dawson integral" begin
    reference = [
        0.0 => 0.0,
        0.1 => 0.09933599239785286,
        0.5 => 0.4244363835020223,
        1.0 => 0.5380795069127684,
        2.0 => 0.30134038892379197,
        5.0 => 0.10213407442427684,
        10.0 => 0.05025384718759853,
        100.0 => 0.005000250037509379,
    ]
    for (x, expected) in reference
        @test isapprox(dawson(x), expected; rtol = 2.0e-14)
    end
    @test dawson(0.0) == 0.0
    @test signbit(dawson(-0.0))
    @test dawson(-2.0) == -dawson(2.0)
    @test isnan(dawson(NaN))
    @test dawson(Inf) == 0.0
end

@testset "Dawson numerical regions" begin
    @test isapprox(dawson(1.0e-10), 1.0e-10; rtol = 2.0e-14)
    @test isapprox(
        FewSpecialFunctions._dawson_small_cf(0.1),
        0.09933599239785286;
        rtol = 2.0e-14,
    )
    @test isapprox(
        FewSpecialFunctions._dawson_middle_series(2.0),
        0.30134038892379197;
        rtol = 2.0e-14,
    )
    @test isapprox(dawson(20.0), 0.02503136792640367; rtol = 2.0e-14)
    @test isapprox(
        FewSpecialFunctions._dawson_large_cf(20.0),
        0.02503136792640367;
        rtol = 2.0e-14,
    )
    @test dawson(1.0e12) == 5.0e-13

    for x in (0.1, 1.0, 10.0)
        h = 1.0e-6
        derivative = (dawson(x + h) - dawson(x - h)) / (2h)
        @test isapprox(derivative, muladd(-2x, dawson(x), 1.0); rtol = 1.0e-7)
    end
end

@testset "Dawson type handling" begin
    @test dawson(1.0f0) isa Float32
    @test dawson(1.0) isa Float64
    @test dawson(big"1.0") isa BigFloat
    @test dawson(1) isa Float64
    @test dawson(1 // 2) isa Float64
end
