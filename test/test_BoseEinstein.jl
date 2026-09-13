using SpecialFunctions: gamma, zeta
using QuadGK: quadgk
using DelimitedFiles: readdlm

@testset "Bose–Einstein integrals" begin
    @test BoseEinsteinIntegralNorm(0.5, -1.0) ≈ 0.4284407345998379 rtol = 3.0e-15
    @test BoseEinsteinIntegralNorm(-4.5, -0.75) ≈ 42.4549731123344429 rtol = 3.0e-15
    @test BoseEinsteinIntegral(1.5, -1.0) ≈ gamma(2.5) * 0.3957280103803376 rtol = 3.0e-15
    @test BoseEinsteinIntegralNorm(1, 0) ≈ pi^2 / 6
    @test BoseEinsteinIntegralNorm(0, -1.0e-20) ≈ -log(1.0e-20)
    @test BoseEinsteinIntegralNorm(-1, -1.0e-20) ≈ 1.0e20
    @test BoseEinsteinIntegralNorm(1, -Inf) == 0.0
    @test BoseEinsteinIntegral(200, -Inf) == 0.0
    @test BoseEinsteinIntegralNorm(-0.5, 0) == Inf
    @test BoseEinsteinIntegral(0, 0) == Inf
    @test BoseEinsteinIntegralNorm(0.5f0, -1.0f0) isa Float32
    @test BoseEinsteinIntegralNorm(1, -1) isa Float64
    @test BoseEinsteinIntegralNorm(big"0.5", big"-1") isa BigFloat
    @test_throws DomainError BoseEinsteinIntegralNorm(0.25, -1)
    @test_throws DomainError BoseEinsteinIntegralNorm(-5, -1)
    @test_throws DomainError BoseEinsteinIntegralNorm(1, 0.1)
    @test_throws DomainError BoseEinsteinIntegralNorm(1, NaN)
    @test_throws DomainError BoseEinsteinIntegralNorm(Inf, -1)
    @test_throws DomainError BoseEinsteinIntegral(-1, -1)

    @testset "Fukushima Tables A.49–A.50" begin
        samples = readdlm(joinpath(@__DIR__, "data", "bose_einstein_fukushima.tsv"), '\t', Float64; skipstart = 1)
        for row in eachrow(samples), (column, η) in enumerate((-1.5, -1.0, -0.75))
            @test BoseEinsteinIntegralNorm(row[1] / 2, η) ≈ row[column + 1] rtol = 3.0e-15
        end
    end

    @testset "Defining integral" begin
        # t = u² removes the integrable origin singularity at k = -1/2.
        for k in (-0.5, 0.0, 0.5, 1.0, 1.5, 5.5, 6.0, 19.5, 20.0), η in (-3.0, -0.1)
            reference, _ = quadgk(0.0, Inf; rtol = 2.0e-13) do u
                t = u^2
                return 2exp((2k + 1) * log(u) - t + η) / (-expm1(η - t))
            end
            @test BoseEinsteinIntegral(k, η) ≈ reference rtol = 5.0e-13
            @test BoseEinsteinIntegralNorm(k, η) ≈ reference / gamma(k + 1) rtol = 5.0e-13
        end
    end

    @testset "Independent high-precision values and switch boundaries" begin
        # Li_{k+1}(exp(η)), summed with Python decimal at 110-digit precision
        # and cross-checked with mpmath 1.3.0 at 130 digits. Includes all 17
        # cutoffs and their adjacent Float64 values, plus both series regions.
        # Near-zero rows (η = -1e-20, -1e-6) use mpmath at 160 digits.
        references = readdlm(joinpath(@__DIR__, "data", "bose_einstein_reference.tsv"), '\t', String; skipstart = 1)
        for row in eachrow(references)
            k, η, reference = parse.(Float64, row)
            k < -4.5 && continue # Lower orders are used by differentiation only.
            @test BoseEinsteinIntegralNorm(k, η) ≈ reference rtol = 3.0e-15
            @test BoseEinsteinIntegralNorm(Float32(k), Float32(η)) ≈ Float32(reference) rtol = 8eps(Float32)
        end
        for bits in (128, 256)
            setprecision(BigFloat, bits) do
                for row in eachrow(references)
                    k, η, reference = parse.(BigFloat, row)
                    k < -4.5 && continue
                    @test BoseEinsteinIntegralNorm(k, η) ≈ reference rtol = 64eps(BigFloat)
                    if k > -1
                        @test BoseEinsteinIntegral(k, η) ≈ gamma(k + 1) * reference rtol = 64eps(BigFloat)
                    end
                end
            end
        end
    end

    @testset "Endpoints, limiting behavior, and large orders" begin
        for k in -4.5:0.5:30
            @test BoseEinsteinIntegralNorm(k, -Inf) == 0
            @test BoseEinsteinIntegralNorm(k, 0) == (k > 0 ? zeta(k + 1) : Inf)
            @test BoseEinsteinIntegralNorm(k, -700) ≈ exp(-700) rtol = 3.0e-15
            @test BoseEinsteinIntegralNorm(k, -2) < BoseEinsteinIntegralNorm(k, -1)
            @test BoseEinsteinIntegralNorm(k, -1.0e-20) > 0
            if k < 0
                @test BoseEinsteinIntegralNorm(k, -1.0e-100) ≈ gamma(-k) * 10.0^(-100k) rtol = 3.0e-15
            elseif k > 0
                @test BoseEinsteinIntegralNorm(k, -1.0e-30) ≈ zeta(k + 1) rtol = 5.0e-15
            end
        end
        for k in (20, 20.5, 32, 53, 100, 200, 1.0e20)
            # At η = -1, summing the first 40 positive terms is an independent
            # reference even for the lowest order here (tail < 1e-50).
            reference = Float64(sum(exp(-big(n)) / big(n)^(big(k) + 1) for n in 1:40))
            @test BoseEinsteinIntegralNorm(k, -1) ≈ reference rtol = 3.0e-15
            @test BoseEinsteinIntegralNorm(k, -1.0e-10) ≈ Float64(BoseEinsteinIntegralNorm(big(k), big"-1e-10")) rtol = 3.0e-15
        end
        for (k, η) in ((20, -750), (170, -745), (200, -1000))
            reference = Float64(gamma(big(k + 1)) * exp(big(η)))
            @test BoseEinsteinIntegral(k, η) ≈ reference rtol = 3.0e-13
        end
        @test BoseEinsteinIntegral(200, -1) == Inf
    end

    @testset "Types, broadcasting, and domain" begin
        for f in (BoseEinsteinIntegral, BoseEinsteinIntegralNorm)
            @test f(Float16(0.5), Float16(-1)) isa Float16
            @test f(0.5f0, -1.0f0) isa Float32
            @test f(1, -1.0f0) isa Float32
            @test f(1 // 2, -1) isa Float64
            @test f(0.5, big"-1") isa BigFloat
            @test f.([0.5, 1.5], [-2.0 -1.0 0.0]) == [f(k, η) for k in (0.5, 1.5), η in (-2.0, -1.0, 0.0)]
            for k in (-Inf, Inf, NaN, 0.1, 1.25), η in (-1.0, 0.0, -Inf)
                @test_throws DomainError f(k, η)
            end
            for η in (NaN, Inf, 1.0e-20)
                @test_throws DomainError f(0.5, η)
            end
        end
        @test @inferred(BoseEinsteinIntegralNorm(1.5, -0.5)) isa Float64
        @test @inferred(BoseEinsteinIntegral(1.5f0, -0.5f0)) isa Float32
    end
end
