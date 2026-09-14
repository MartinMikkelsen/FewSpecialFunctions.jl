include("data/whittaker.jl")

@testset "Whittaker functions" begin
    for (κ, μ, z, m, w, dm, dw) in WHITTAKER_REFERENCES
        @test WhittakerM(κ, μ, z) ≈ m rtol = 2.0e-12
        @test WhittakerW(κ, μ, z) ≈ w rtol = 2.0e-12
        @test dWhittakerM(κ, μ, z) ≈ dm rtol = 2.0e-12
        @test dWhittakerW(κ, μ, z) ≈ dw rtol = 2.0e-12
    end

    # Elementary cases catch normalization and parameter-shift errors.
    for z in (0.1, 1.0, 100.0, 1 + 2im)
        @test WhittakerM(0, 0.5, z) ≈ 2sinh(z / 2) rtol = 2.0e-13
        @test WhittakerW(0, 0.5, z) ≈ exp(-z / 2) rtol = 2.0e-13
        @test dWhittakerW(0, 0.5, z) ≈ -exp(-z / 2) / 2 rtol = 2.0e-13
        @test WhittakerW(1, 0.5, z) ≈ z * exp(-z / 2) rtol = 2.0e-13
        @test WhittakerW(0.2, 0.3, z) ≈ WhittakerW(0.2, -0.3, z) rtol = 2.0e-13
    end
    for f in (WhittakerM, WhittakerW, dWhittakerM, dWhittakerW)
        @test f(0, 1, 2) isa Float64
        @test f(0.2f0, 0.3f0, 1.0f0) isa Float32
        @test f(0.2f0, 0.3f0, 1.0f0) ≈ Float32(f(Float64(0.2f0), Float64(0.3f0), 1.0)) rtol = 2eps(Float32)
        @test f(0.2, 0.3, 1 + 2im) isa ComplexF64
        @test_throws DomainError f(0.2, 0.3, 0.0)
        @test_throws DomainError f(0.2, 0.3, -1.0)
        @test_throws DomainError f(Inf, 0.3, 1.0)
        @test_throws DomainError f(0.2, NaN, 1.0)
        @test_throws DomainError f(0.2, 0.3, Inf)
    end
    @test_throws DomainError WhittakerM(0.2, -0.5, 1.0)
    @test_throws DomainError dWhittakerM(0.2, -1.0, 1.0)
    @test WhittakerW(0, 0.5, 1.0e6) == 0.0
    @test dWhittakerW(0, 0.5, 1.0e6) == -0.0
    for f in (WhittakerM, WhittakerW)
        @test f(0.2, 0.3, complex(-2.0, -0.0)) ≈ conj(f(0.2, 0.3, complex(-2.0, 0.0))) rtol = 2.0e-13
    end
    setprecision(192) do
        # Independent mpmath values at 90 decimal digits for decimal inputs.
        κ, μ, z = parse.(BigFloat, ("0.7", "0.4", "1.3"))
        @test WhittakerM(κ, μ, z) ≈ parse(BigFloat, "0.79162087513607241909724796869043879582406712375568991994834961183166660989188673526736526") rtol = big"1e-50"
        @test WhittakerW(κ, μ, z) ≈ parse(BigFloat, "0.678444373413397358673290850662677327534407158149638528762132925320949920731236328707314511") rtol = big"1e-50"
        @test dWhittakerW(κ, μ, z) ≈ parse(BigFloat, "-0.00966302157752613773418353922599935922958468069562053600458587133605871565288565160193487046") rtol = big"1e-50"
        @test WhittakerW(κ, μ, z + im) isa Complex{BigFloat}
        @test WhittakerM(big(0), big"0.5", big(2)) ≈ 2sinh(big(1)) rtol = big"1e-50"
        @test WhittakerW(big(0), big"0.5", big(2)) ≈ exp(-big(1)) rtol = big"1e-50"
        @test precision(WhittakerW(big"0.2", big"0.3", big(1))) == 192
        @test precision(BigFloat) == 192
    end
    # Hitting the work limit must reject an incomplete asymptotic sum.
    setprecision(8192) do
        @test FewSpecialFunctions._whittaker_w_asymptotic(BigFloat(0), BigFloat(0), BigFloat(10000)) === nothing
    end
end
