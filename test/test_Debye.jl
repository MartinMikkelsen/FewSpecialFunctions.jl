using DelimitedFiles
using QuadGK: quadgk
using SpecialFunctions: zeta
using Test

@testset "Debye" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "debye_test.txt"))

    @testset "Debye_function vs DebyeFunctions.jl" begin
        for r in 1:size(data, 1)
            x, n, Q_ref = data[r, :]
            @test FewSpecialFunctions.debye_function(n, 1.0, x) ≈ Q_ref atol = 1.0e-4
        end
    end

    @test FewSpecialFunctions.debye_function(2.0, 1.0, 5.0) ≈ 0.172329034857624782145 atol = 1.0e-6
    @test FewSpecialFunctions.debye_function(5.0, 1.0, 4.5) ≈ 0.10164118339698890968 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(7.0, 1.0, 0.8) ≈ 0.69112406526865230673 atol = 1.0e-14
    @test FewSpecialFunctions.debye_function(9.0, 1.0, 3.4) ≈ 0.15413773867789254146 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(12.0, 1.0, 5.4) ≈ 0.03618849233828133 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(15.0, 1.0, 2.4) ≈ 0.2661409156647294951955 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(20.0, 1.0, 1.24) ≈ 0.523361585088859680745 atol = 1.0e-11
    @test FewSpecialFunctions.debye_function(25.0, 1.0, 4.2) ≈ 0.07296496706218587 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(30.0, 1.0, 3.42) ≈ 0.1258426106590655660781 atol = 1.0e-13

end

@testset "Debye function method variants" begin
    # Test real number conversion
    @test FewSpecialFunctions.debye_function(2, 1, 5) ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 5.0)
    @test FewSpecialFunctions.debye_function(Int8(3), Float32(1.0), 2) ≈ FewSpecialFunctions.debye_function(3.0, 1.0, 2.0)

    # Test array input for x
    x_array = [1.0, 2.0, 3.0]
    result = FewSpecialFunctions.debye_function(2.0, 1.0, x_array)
    @test length(result) == length(x_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 1.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 2.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)

    # Test array input for β
    β_array = [1.0, 1.5, 2.0]
    result = FewSpecialFunctions.debye_function(2.0, β_array, 3.0)
    @test length(result) == length(β_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.5, 3.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(2.0, 2.0, 3.0)

    # Test array input for n
    n_array = [1.0, 2.0, 3.0]
    result = FewSpecialFunctions.debye_function(n_array, 1.0, 3.0)
    @test length(result) == length(n_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(1.0, 1.0, 3.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(3.0, 1.0, 3.0)

    # Test default n=1.0
    @test FewSpecialFunctions.debye_function(1.5, 2.0) ≈ FewSpecialFunctions.debye_function(1.0, 1.5, 2.0)
end

@testset "debye_function validation" begin
    @test_throws ArgumentError FewSpecialFunctions.debye_function(1.0, 1.0, -1.0)
    @test FewSpecialFunctions.debye_function(1.0, 1.0, 0.0) == 1.0  # limit at x=0 is 1
    @test_throws ArgumentError FewSpecialFunctions.debye_function(1.0, -1.0, 1.0)
    @test_throws ArgumentError FewSpecialFunctions.debye_function(-1.0, 1.0, 1.0)
end

@testset "Debye convergence and generalized limits" begin
    for args in ((Float16(2), Float16(1), Float16(1)), (Float16(2), Int8(1), Float16(1)))
        result = debye_function(args...)
        @test result isa Float16
        @test result ≈ 0.7078784756278294 rtol = 8eps(Float16)
        @test_throws ArgumentError debye_function(args...; tol = 0)
    end

    # Integrating the Bernoulli expansion gives D₂(x) = 1 - x/3 + x²/24 - x⁴/2160 + O(x⁶).
    for T in (Float32, Float64, BigFloat), x in (T(1.0e-6), T(0.001))
        expected = one(T) - x / 3 + x^2 / 24 - x^4 / 2160
        result = debye_function(T(2), one(T), x)
        @test result isa T
        @test result ≈ expected rtol = max(T(1.0e-20), 16eps(T))
    end

    # Independent high-precision integration in t = v², without the production transform.
    setprecision(BigFloat, 192) do
        for (n, β, x) in ((2, 1, 1), (2, 2, 0.01), (1, 1.5, 2), (3, 0.5, 4))
            nb, βb, xb = BigFloat(n), BigFloat(β), BigFloat(x)
            integral, error = quadgk(
                v -> 2v^(2nb + 1) / expm1(v^2)^βb,
                zero(xb), sqrt(xb); rtol = big"1e-45"
            )
            reference = nb * integral / xb^nb
            @test error <= big"1e-45" * integral
            @test debye_function(Float64(n), Float64(β), Float64(x)) ≈ reference rtol = 2.0e-14
            @test debye_function(nb, βb, xb) ≈ reference rtol = big"1e-34"
        end
    end

    # At large x the omitted upper tails are exponentially small.
    @test debye_function(2.0, 1.0, 1.0e6) ≈ 4zeta(3.0) / 1.0e12 rtol = 2.0e-14
    @test debye_function(2.0, 1.0, 1.0e150) ≈ 4zeta(3.0) / 1.0e300 rtol = 2.0e-14
    @test debye_function(2.0, 2.0, 1.0e6) ≈ 4(zeta(2.0) - zeta(3.0)) / 1.0e12 rtol = 2.0e-14

    @test debye_function(2.0, 0.5, 0.0) == 0.0
    @test debye_function(2.0, 1.0, 0.0) == 1.0
    @test debye_function(2.0, 2.0, 0.0) == Inf
    @test debye_function(1.0e-20, 1.0, 1.0) ≈ 1.0
    @test debye_function(1.0e20, 1.0e20, 0.0) == Inf
    @test debye_function(2.0, 1.0, Inf) == 0.0
    for β in (3.0, 4.0), x in (0.0, 1.0)
        @test_throws ArgumentError debye_function(2.0, β, x)
    end
    for tol in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError debye_function(2.0, 1.0, 1.0; tol)
    end
    for max_terms in (0, -1, 1.5)
        @test_throws ArgumentError debye_function(2.0, 1.0, 1.0; max_terms)
    end
    @test_throws ErrorException debye_function(2.0, 1.0, 1.0; max_terms = 1)
    @test_throws ErrorException debye_function(1.0e20, 1.0, 1.0)
    @test debye_function(2.0, 1.0, 1.0; tol = 1.0e-8, max_terms = 100) ≈ 0.7078784756278294 rtol = 1.0e-8
end
