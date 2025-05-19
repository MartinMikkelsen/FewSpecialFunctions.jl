using SpecialFunctions

import FewSpecialFunctions: Q, c_μ, lnA, half_ζ2 ,ζ ,theta_over_sin,theta_prime_sin,ρ,r,r_prime_sin,f,ψ, MarcumQ_small_x, MarcumQ_large_xy, MarcumQ_recurrence, MarcumQ_large_M, MarcumQ_quadrature, MarcumQ_modified, f1_f2

@testset "Marcum Q-function" begin

    # assume marcumq_test.txt lives in test/data/
    data = open(readdlm, joinpath(@__DIR__, "data", "marcumq_test.txt"))

    @testset "MarcumQ vs MATLAB reference" begin
        @test size(data, 2) == 4

        for r in 1:size(data, 1)
            M, a, b, Q_ref = data[r, :]

            # test the standard MarcumQ(M,a,b)
            @test FewSpecialFunctions.MarcumQ(M, a, b) ≈ Q_ref rtol = 1e-10
        end
    end
    

    # Acta Univ. Sapientiae, Mathematica, 3, 1 (2011) 60–7
    @test FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) ≈ 0.838249985438908 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 0.2, 0.6) ≈ 0.999998670306184 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 0.2, 0.6) ≈  0.999999999927717 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 1.2, 1.6) ≈ 0.501536568390858 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 1.2, 1.6) ≈ 0.994346394491553 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 1.2, 1.6) ≈   0.999944937223540 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9

    # Compared to marcumq (https://se.mathworks.com/help/signal/ref/marcumq.html)

    @test FewSpecialFunctions.MarcumQ(5, 0.47, 4.85) ≈ 0.0106766402997493 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 1.46, 4.00) ≈ 0.211798804811782  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 1.27, 4.58) ≈ 0.000931257801666407 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 2.88, 3.28) ≈ 0.773155207859263 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 2.55, 4.67) ≈  0.024112315799424 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 2.27, 3.72) ≈ 0.400088995953665 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(2, 1.97, 0.86) ≈ 0.990345203236692 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 0.10, 1.38) ≈ 0.983869651076909 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 0.29, 4.12) ≈ 0.00028475422733874 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 0.95, 4.75) ≈ 0.00899422673877906 atol = 1e-9

end
@testset "MarcumQ overloads" begin
    # Test for array inputs
    @testset "Array inputs" begin
        # Test when b is an array
        b_arr = [0.6, 1.6, 2.6]
        q_b_arr = FewSpecialFunctions.MarcumQ(1, 0.2, b_arr)
        @test q_b_arr isa Vector
        @test length(q_b_arr) == length(b_arr)
        @test q_b_arr[1] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test q_b_arr[2] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 1.6) atol = 1e-9
        @test q_b_arr[3] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 2.6) atol = 1e-9
        
        # Test when a is an array
        a_arr = [0.2, 1.2, 2.2]
        q_a_arr = FewSpecialFunctions.MarcumQ(1, a_arr, 0.6)
        @test q_a_arr isa Vector
        @test length(q_a_arr) == length(a_arr)
        @test q_a_arr[1] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test q_a_arr[2] ≈ FewSpecialFunctions.MarcumQ(1, 1.2, 0.6) atol = 1e-9
        @test q_a_arr[3] ≈ FewSpecialFunctions.MarcumQ(1, 2.2, 0.6) atol = 1e-9
        
        # Test when μ is an array
        μ_arr = [1.0, 5.0, 7.7]
        q_μ_arr = FewSpecialFunctions.MarcumQ(μ_arr, 0.2, 0.6)
        @test q_μ_arr isa Vector
        @test length(q_μ_arr) == length(μ_arr)
        @test q_μ_arr[1] ≈ FewSpecialFunctions.MarcumQ(1.0, 0.2, 0.6) atol = 1e-9
        @test q_μ_arr[2] ≈ FewSpecialFunctions.MarcumQ(5.0, 0.2, 0.6) atol = 1e-9
        @test q_μ_arr[3] ≈ FewSpecialFunctions.MarcumQ(7.7, 0.2, 0.6) atol = 1e-9
    end
    
    # Test for the default μ=1 case
    @testset "Default μ=1" begin
        @test FewSpecialFunctions.MarcumQ(0.2, 0.6) ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test FewSpecialFunctions.MarcumQ(1.2, 1.6) ≈ FewSpecialFunctions.MarcumQ(1, 1.2, 1.6) atol = 1e-9
        @test FewSpecialFunctions.MarcumQ(2.2, 2.6) ≈ FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) atol = 1e-9
    end
end

@testset "dQdb - Marcum Q derivative" begin
    @test FewSpecialFunctions.dQdb(1, 1.0, 2.0) < 0  
    
    function numerical_dQdb(M, a, b, h=1e-6)
        return (FewSpecialFunctions.MarcumQ(M, a, b + h) - FewSpecialFunctions.MarcumQ(M, a, b - h)) / (2h)
    end
    
    test_cases = [
        (1, 1.0, 2.0), 
        (2, 1.5, 2.5), 
        (3, 0.5, 1.5),
        (4, 2.0, 3.0)
    ]
    
    for (M, a, b) in test_cases
        analytical = FewSpecialFunctions.dQdb(M, a, b)
        numerical = numerical_dQdb(M, a, b)
        @test analytical ≈ numerical rtol = 1e-4
    end
    
    # Edge cases
    @test FewSpecialFunctions.dQdb(1, 0.1, 0.1) ≈ numerical_dQdb(1, 0.1, 0.1) rtol = 1e-4  # Small values
    @test FewSpecialFunctions.dQdb(10, 5.0, 5.0) ≈ numerical_dQdb(10, 5.0, 5.0) rtol = 1e-3  # Larger order
    
    # Error conditions
    @test_throws AssertionError FewSpecialFunctions.dQdb(0, 1.0, 2.0)  # M must be ≥ 1
end

@testset "Helper functions for MarcumQ" begin
    @testset "Q function (regularized incomplete gamma)" begin
        # Test against known values
        @test FewSpecialFunctions.Q(1.0, 0.0) ≈ 1.0 atol = 1e-14
        @test FewSpecialFunctions.Q(1.0, 1.0) ≈ exp(-1.0) atol = 1e-14
        @test FewSpecialFunctions.Q(2.0, 2.0) ≈ exp(-2.0) * (1.0 + 2.0) atol = 1e-14
        @test FewSpecialFunctions.Q(3.0, 3.0) ≈ exp(-3.0) * (1.0 + 3.0 + 9.0/2.0) atol = 1e-14
        
        # Compare with direct computation for large values
        μ, x = 5.0, 10.0
        @test FewSpecialFunctions.Q(μ, x) ≈ gamma_inc(μ, x)[2] atol = 1e-14
    end

    @testset "c_μ function (Bessel I ratio)" begin
        # Test specific values
        for μ in [0.5, 1.0, 1.5, 2.0]
            for ξ in [0.1, 0.5, 1.0, 2.0, 5.0]
                direct_ratio = besseli(μ, ξ) / besseli(μ-1, ξ)
                @test FewSpecialFunctions.c_μ(μ, ξ) ≈ direct_ratio atol = 1e-14
            end
        end
        
        # Test limit behavior for large ξ
        # For large ξ, c_μ(μ, ξ) ≈ 1 - (μ-0.5)/ξ
        large_ξ = 100.0
        for μ in [1.0, 2.0, 3.0]
            approx = 1.0 - (μ-0.5)/large_ξ
            @test FewSpecialFunctions.c_μ(μ, large_ξ) ≈ approx rtol = 1e-3
        end
    end
end

using SpecialFunctions: besseli

@testset "FewSpecialFunctions dQdb Tests" begin

    # Scalar input tests
    @testset "dQdb Scalar Tests" begin
        M, a, b = 1, 2.0, 3.0
        coeff = b^M / a^(M-1)
        expected = -coeff * exp(-(a^2 + b^2) / 2) * besseli(M-1, a * b)
        result = dQdb(M, a, b)
        @test isapprox(result, expected, atol=1e-6)

        # Convenience method: M=1
        a, b = 2.0, 3.0
        result = dQdb(a, b)
        @test isapprox(result, dQdb(1, a, b), atol=1e-6)
    end

    # Array input tests
    @testset "dQdb Array Tests" begin

        # Array handling for b parameter
        M, a, b = 1, 2.0, [0.0, 3.0]
        results = dQdb(M, a, b)
        expected1 = -exp(-4.0)
        expected2 = -3.0 * exp(-6.5) * besseli(0, 6.0)
        @test isapprox(results[2], expected2, atol=1e-6)

        # Array handling for a parameter
        M, a, b = 1, [1.0, 2.0], 3.0
        results = dQdb(M, a, b)
        expected1 = -3.0 * exp(-5.0) * besseli(0, 3.0)
        expected2 = -3.0 * exp(-6.5) * besseli(0, 6.0)
        @test isapprox(results[1], expected1, atol=1e-6)
        @test isapprox(results[2], expected2, atol=1e-6)

        # Array handling for M parameter
        M, a, b = [1, 2], 2.0, 3.0
        results = dQdb(M, a, b)
        expected1 = -3.0 * exp(-6.5) * besseli(0, 6.0)
        expected2 = -9.0 / 2.0 * exp(-6.5) * besseli(1, 6.0)
        @test isapprox(results[1], expected1, atol=1e-6)
        @test isapprox(results[2], expected2, atol=1e-6)
    end

end

@testset "Q (upper incomplete gamma)" begin
    # Γ(μ,x) for μ=1 is e^{-x}
    @test Q(1.0, 2.0) ≈ exp(-2.0) atol=1e-9
    # Γ(2,x) = e^{-x}(x+1)
    @test Q(2.0, 3.0) ≈ exp(-3.0)*(3.0 + 1.0) atol=1e-9
    # at x=0, Q(μ,0) == Γ(μ)
    μ = 3.5
    @test Q(μ, 0.0) ≈ 1.0 atol=1e-9
    @test Q(μ, 100.0) ≈ 0.0 atol=1e-9
end

@testset "c_μ (continued‐fraction ratio for scaled I_μ)" begin
    for μ in (0.0, 0.5, 1.0, 2.5)
        for ξ in (0.1, 1.0, 5.0)
            @test c_μ(μ, ξ) ≈ besseli(μ, ξ)/besseli(μ-1, ξ) atol=1e-9
        end
    end
end

@testset "lnA (log Aₙ from Eq. 32)" begin
    # n=0 should give lnA=0
    @test lnA(0, 1.3) ≈ 0.0 atol=1e-9
    # for n=1, μ=1: A₁ = 3/8
    val = lnA(1, 1.0)
    @test val ≈ log(3/8) atol=1e-9
    # for n=2, μ=0.5: compute expected directly
    μ, n = 0.5, 2
    expected = loggamma(μ+0.5+n) - loggamma(μ+0.5-n) - n*log(2.0) - loggamma(n+1)
    @test lnA(n, μ) ≈ expected atol=1e-9
end

@testset "half_ζ2 (ζ²/2 from Eq. 84)" begin
    # small δ branch: δ = y - x - 1 ≈ 0 ⇒ half_ζ2 → 0
    for x in (0.5, 1.0, 2.0)
        y = x + 1.0
        @test half_ζ2(x, y) ≈ 0.0 atol=1e-9
    end
    # general branch
    x, y = 1.2, 3.4
    δ = y - x - 1.0
    @test abs(δ) > 1e-3
    r = sqrt(1 + 4*x*y)
    expected = x + y - r + log((1 + r)/(2*y))
    @test half_ζ2(x, y) ≈ expected atol=1e-9
end


@testset "ζ (Eq. 56)" begin
    @test ζ(0.5, 1.5) ≈ copysign(sqrt(2 * half_ζ2(0.5, 1.5)), 0.5 + 1.0 - 1.5) atol=1e-9
end

@testset "θ/sin(θ) and θ′sin(θ) for small θ" begin
    @test theta_over_sin(0.0) ≈ 1.0 atol=1e-9
    @test theta_over_sin(1e-6) ≈ 1.0 + (1e-6)^2/6 atol=1e-9
    @test theta_over_sin(0.1) ≈ 0.1/sin(0.1) atol=1e-9
    @test theta_prime_sin(0.0) ≈ 0.0 atol=1e-9
    @test theta_prime_sin(1e-6) ≈ (1e-6)^2/3 atol=1e-9
    @test theta_prime_sin(0.1) ≈ 1.0 - 0.1/tan(0.1) atol=1e-9
end

@testset "ρ from eq. (98)" begin
    @test ρ(0.0, 1.0) ≈ 1.0 atol=1e-9
    @test ρ(1.0, 0.0) ≈ 1.0 atol=1e-9
    @test ρ(1.0, 1.0) ≈ sqrt(2.0) atol=1e-9
end

@testset "r and r′sin(θ) from eq. (96)" begin
    θ, y, ξ = 0.1, 2.0, 0.5
    tos = theta_over_sin(θ)
    xs = ξ / tos
    expected_r = (1.0 + sqrt(1.0 + xs^2)) * tos / (2.0 * y)
    @test r(θ, y, ξ) ≈ expected_r atol=1e-9
    expected_r_prime = (1.0 + 1.0 / sqrt(1.0 + xs^2)) * theta_prime_sin(θ) / (2.0 * y)
    @test r_prime_sin(θ, y, ξ) ≈ expected_r_prime atol=1e-9
end

@testset "f from eq. (96)" begin
    θ, y, ξ = 0.1, 2.0, 0.5
    r0 = r(θ, y, ξ)
    d = r0 - cos(θ)
    expected_f = (r_prime_sin(θ, y, ξ) - d * r0) / (d^2 + sin(θ)^2)
    @test f(θ, y, ξ) ≈ expected_f atol=1e-9
end

@testset "ψ from eq. (97)" begin
    θ, ξ = 0.1, 0.5
    tos = theta_over_sin(θ)
    rv = ρ(tos, ξ)
    expected_ψ = cos(θ) * rv - sqrt(1.0 + ξ^2) - log((tos + rv) / (1.0 + sqrt(1.0 + ξ^2)))
    @test ψ(θ, ξ) ≈ expected_ψ atol=1e-9
end

@testset "f1_f2 (boundaries for recurrence relation)" begin
    # Test specific values
    x, M = 1.0, 2.0
    f1, f2 = f1_f2(x, M)
    @test f1 ≈ x + M - sqrt(4*x + 2*M) atol=1e-9
    @test f2 ≈ x + M + sqrt(4*x + 2*M) atol=1e-9
    
    # Test multiple values
    test_cases = [
        (0.5, 1.0),
        (2.0, 3.0),
        (10.0, 5.0),
        (0.1, 0.2)
    ]
    
    for (x, M) in test_cases
        f1, f2 = f1_f2(x, M)
        expected_f1 = x + M - sqrt(4*x + 2*M)
        expected_f2 = x + M + sqrt(4*x + 2*M)
        @test f1 ≈ expected_f1 atol=1e-9
        @test f2 ≈ expected_f2 atol=1e-9
        
        # Verify that f1 < f2
        @test f1 < f2
        
        # Verify that f1, f2 represent a valid interval
        @test f2 - f1 ≈ 2*sqrt(4*x + 2*M) atol=1e-9
    end
    
    # Test with different numeric types
    x, M = 1.0f0, 2.0f0  # Float32
    f1, f2 = f1_f2(x, M)
    @test typeof(f1) == Float32
    @test typeof(f2) == Float32
    @test f1 ≈ x + M - sqrt(4*x + 2*M) atol=1e-6
    @test f2 ≈ x + M + sqrt(4*x + 2*M) atol=1e-6
end

