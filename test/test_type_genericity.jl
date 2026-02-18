using Test
using SpecialFunctions

@testset "Type genericity" begin

    @testset "Debye Float32" begin
        result = FewSpecialFunctions.debye_function(2.0f0, 1.0f0, 5.0f0)
        @test result isa Float32
        @test isapprox(result, 0.172329f0; atol = 1.0f-3)

        # Integer inputs promote to Float64
        result_int = FewSpecialFunctions.debye_function(2, 1, 5)
        @test result_int isa Float64

        # Mixed types promote
        result_mix = FewSpecialFunctions.debye_function(2.0f0, 1.0, 5.0f0)
        @test result_mix isa Float64
    end

    @testset "Debye BigFloat" begin
        result = FewSpecialFunctions.debye_function(big"2.0", big"1.0", big"5.0")
        @test result isa BigFloat
        @test isapprox(result, big"0.172329034857624782145"; atol = big"1e-5")
    end

    @testset "FermiDirac Float32" begin
        result = FewSpecialFunctions.FermiDiracIntegral(0.0f0, 0.0f0)
        @test result isa Float32
        @test isapprox(result, 0.6931472f0; atol = 1.0f-4)

        # j=0 branch: log(1 + exp(x))
        result2 = FewSpecialFunctions.FermiDiracIntegral(0.0f0, 1.0f0)
        @test result2 isa Float32
        @test isapprox(result2, Float32(log(1 + exp(1.0))); atol = 1.0f-4)
    end

    @testset "FermiDirac generic approximation branch" begin
        # The generic approximation (j not a half-integer) should work with Float32
        result = FewSpecialFunctions.FermiDiracIntegral(1.0f0, 1.0f0)
        @test result isa Float32
        @test isfinite(result)
    end

    @testset "Parabolic cylinder Float32" begin
        result_U = FewSpecialFunctions.U(0.0f0, 1.0f0)
        @test result_U isa Float32
        @test isapprox(result_U, Float32(FewSpecialFunctions.U(0.0, 1.0)); atol = 1.0f-4)

        result_V = FewSpecialFunctions.V(0.0f0, 1.0f0)
        @test result_V isa Float32
        @test isapprox(result_V, Float32(FewSpecialFunctions.V(0.0, 1.0)); atol = 1.0f-4)

        result_W = FewSpecialFunctions.W(0.0f0, 1.0f0)
        @test result_W isa Float32
        @test isapprox(result_W, Float32(FewSpecialFunctions.W(0.0, 1.0)); atol = 1.0f-4)
    end

    @testset "Parabolic cylinder derivatives Float32" begin
        result_dU = FewSpecialFunctions.dU(0.0f0, 1.0f0)
        @test result_dU isa Float32
        @test isapprox(result_dU, Float32(FewSpecialFunctions.dU(0.0, 1.0)); atol = 1.0f-4)

        result_dV = FewSpecialFunctions.dV(0.0f0, 1.0f0)
        @test result_dV isa Float32
        @test isapprox(result_dV, Float32(FewSpecialFunctions.dV(0.0, 1.0)); atol = 1.0f-4)

        result_dW = FewSpecialFunctions.dW(0.0f0, 1.0f0)
        @test result_dW isa Float32
        @test isapprox(result_dW, Float32(FewSpecialFunctions.dW(0.0, 1.0)); atol = 1.0f-4)
    end

    @testset "Parabolic cylinder integer promotion" begin
        # Integer inputs should promote to Float64
        @test FewSpecialFunctions.U(0, 1) isa Float64
        @test FewSpecialFunctions.V(0, 1) isa Float64
        @test FewSpecialFunctions.W(0, 1) isa Float64
    end

    @testset "MarcumQ Float32" begin
        import FewSpecialFunctions: f1_f2, MarcumQ_large_M, MarcumQ_recurrence

        x, M = 1.0f0, 2.0f0
        f1, f2 = f1_f2(x, M)
        @test typeof(f1) == Float32
        @test typeof(f2) == Float32

        # MarcumQ with Float32 inputs
        result = FewSpecialFunctions.MarcumQ(1.0f0, 0.2f0, 0.6f0)
        @test result isa Float32
        @test isapprox(result, 0.83825f0; atol = 1.0f-3)
    end

    @testset "MarcumQ integer promotion" begin
        # Integer inputs should be promoted
        result = FewSpecialFunctions.MarcumQ(1, 0.2, 0.6)
        @test result isa Float64
    end

    @testset "Clausen type handling" begin
        # f_n with Float32
        result = FewSpecialFunctions.f_n(2, 1, Float32(π) / 4)
        @test result isa Float32

        # Clausen with Float64 (default)
        result = FewSpecialFunctions.Clausen(1, 0.5)
        @test result isa Float64

        # Integer inputs get promoted
        result = FewSpecialFunctions.Clausen(2, 1)
        @test result isa Float64
    end

    @testset "Fresnel type preservation" begin
        # Float64
        C64, S64, E64 = FewSpecialFunctions.fresnel(1.0)
        @test C64 isa Float64
        @test S64 isa Float64

        # Float32 inputs work and produce correct values
        # (erf from SpecialFunctions may promote to Float64 internally)
        C32, S32, E32 = FewSpecialFunctions.fresnel(1.0f0)
        @test isapprox(C32, C64; atol = 1.0e-5)
        @test isapprox(S32, S64; atol = 1.0e-5)
    end

    @testset "Coulomb type handling" begin
        # These use ::Number dispatch, which is maximally generic
        # Integer inputs should work
        @test FewSpecialFunctions.η(2, 3) isa Real
        @test FewSpecialFunctions.C(0, 0.0) isa Real

        # Float32 inputs
        val = FewSpecialFunctions.F(0.0f0, 0.0f0, 1.0f0)
        @test isfinite(val)
    end

end
