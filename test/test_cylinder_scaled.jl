include("data/cylinder_scaled.jl")

@testset "Parabolic cylinder D and scaled functions" begin
    for x in (-10.0, -1.0, 0.0, 1.0, 10.0)
        @test ParabolicCylinderD(0, x) ≈ exp(-x^2 / 4)
        @test ParabolicCylinderD(1, x) ≈ x * exp(-x^2 / 4)
        @test ParabolicCylinderD(2, x) ≈ (x^2 - 1) * exp(-x^2 / 4) rtol = 2.0e-13 atol = 1.0e-28
        @test dParabolicCylinderD(0, x) ≈ -x * exp(-x^2 / 4) / 2
    end
    for (a, x, u, v) in CYLINDER_SCALED_REFERENCES
        @test U_scaled(a, x) ≈ u rtol = 2.0e-12
        @test V_scaled(a, x) ≈ v rtol = 2.0e-12
        @test ParabolicCylinderD_scaled(-a - 0.5, x) ≈ u rtol = 2.0e-12
    end
    @test U(0.0, 60.0) == 0.0
    @test isinf(V(0.0, 60.0))
    for f in (U_scaled, V_scaled, ParabolicCylinderD_scaled)
        @test f(1, 2) isa Float64
        @test f(1.0f0, 2.0f0) isa Float32
        @test_throws DomainError f(0.0, -1.0)
        @test_throws DomainError f(Inf, 1.0)
        @test_throws DomainError f(0.0, Inf)
        @test_throws DomainError f(NaN, 1.0)
    end
    @test ParabolicCylinderD(UInt(1), 2) ≈ 2exp(-1)
    @test dParabolicCylinderD(UInt(0), 2) ≈ -exp(-1)
    for f in (ParabolicCylinderD, dParabolicCylinderD)
        @test_throws DomainError f(0.0, Inf)
        @test_throws DomainError f(Inf, 1.0)
    end
    setprecision(192) do
        @test U_scaled(big(0), big(60)) ≈ parse(BigFloat, "0.129086005176834271877173624296745546111326791811649956225902028350525618591721614376756922") rtol = big"1e-50"
        @test precision(U_scaled(big(0), big(2))) == 192
        @test precision(BigFloat) == 192
    end
    setprecision(8192) do
        @test FewSpecialFunctions._cylinder_v_asymptotic_scaled(BigFloat(0), BigFloat(100)) === nothing
    end
end
