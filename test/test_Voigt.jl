using FewSpecialFunctions
using Test
using QuadGK

function voigt_quadgk_reference(x::Real, y::Real)
    return setprecision(BigFloat, 256) do
        x_big, y_big = BigFloat(x), BigFloat(y)
        value, = quadgk(
            t -> exp(-t^2 / 4 - y_big * t) * cos(x_big * t),
            BigFloat(0),
            BigFloat(30);
            rtol = BigFloat("1e-60"),
        )
        return value / sqrt(BigFloat(π))
    end
end

@testset "Voigt function" begin
    @test hasmethod(voigt, Tuple{Float64, Float64})

    @test voigt(0.0, 0.0) == 1.0
    @test voigt(2.0, 0.0) == exp(-4.0)
    @test_throws DomainError voigt(1.0, -0.1)
    @test voigt(2.0, 0.5) == voigt(-2.0, 0.5)
    @test voigt(1.0f0, 0.5) isa Float64

    # 256-bit QuadGK evaluation of the exact transformed integral on [0, 30].
    reference = [
        (0.0, 0.1),
        (1.0, 0.1),
        (3.0, 0.01),
        (6.0, 4.0),
    ]
    for (x, y) in reference
        @test isapprox(voigt(x, y), Float64(voigt_quadgk_reference(x, y)); rtol = 2.0e-8)
    end

    narrow_y_reference = [
        (1e-12, 1e-10),
        (π / 12 + 1e-12, 1e-10),
    ]
    for (x, y) in narrow_y_reference
        @test isapprox(voigt(x, y), Float64(voigt_quadgk_reference(x, y)); rtol = 2.0e-12)
    end

    @test isapprox(voigt(1.0, 1e200), inv(sqrt(π) * 1e200); rtol = 2.0e-12)

    @test isnan(voigt(NaN, 0.5))
    @test isnan(voigt(0.5, NaN))
    @test voigt(Inf, 0.5) == 0.0
    @test voigt(0.5, Inf) == 0.0
end
