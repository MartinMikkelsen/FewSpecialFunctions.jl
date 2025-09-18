using HypergeometricFunctions
using SpecialFunctions
@testset "Coulomb" begin

    @testset "CoulombF" begin

        data = open(readdlm, joinpath(@__DIR__, "data", "CoulombF.txt"))

        for r in 1:size(data, 1)
            row = data[r, :]
            ℓ = row[1]
            η = row[2]
            ρ = row[3]
            expected = row[4]

            @test real.(FewSpecialFunctions.F(ℓ, η, ρ)) ≈ expected atol = 1.0e-8

        end

    end

    @testset "CoulombG" begin

        data = open(readdlm, joinpath(@__DIR__, "data", "CoulombG.txt"))

        for r in 1:size(data, 1)
            row = data[r, :]
            ℓ = row[1]
            η = row[2]
            ρ = row[3]
            expected = row[4]

            @test real.(FewSpecialFunctions.G(ℓ, η, ρ)) ≈ expected atol = 1.0e-5

        end

    end
    @testset "CoulombF_imag" begin

        data = open(readdlm, joinpath(@__DIR__, "data", "CoulombF.txt"))

        for r in 1:size(data, 1)
            row = data[r, :]
            ℓ = row[1]
            η = row[2]
            ρ = row[3]
            @test real.(FewSpecialFunctions.F(ℓ, η, ρ)) ≈ real.(FewSpecialFunctions.F_imag(ℓ, η, ρ)) atol = 1.0e-7

        end
    end

    @testset "η" begin
        @test η(2.0, 0.5) == 1 / (2 * 0.5)
        @test η(4.0) ≈ 1 / 2.0 atol = 1.0e-10
    end

    @testset "C normalization" begin
        # for η=0, ℓ=0: C = 1
        @test C(0, 0.0) ≈ 1.0 atol = 1.0e-10
        # for η=0, ℓ=1: C = 2^1 * |Γ(2)/Γ(4)| = 2 * (1!)/(3!) = 2/6 = 1/3
        @test C(1, 0.0) ≈ 1 // 3 atol = 1.0e-10
    end

    @testset "θ phase" begin
        # θ is real for real inputs
        val = θ(1, 0.2, 1.5)
        @test isa(val, Real)
    end

    @testset "F regular Coulomb function" begin
        # for η=0, ℓ=0: F(0,0,ρ) = sin(ρ)
        for ρ in (0.1, 1.0, 2.5, π / 3)
            @test real(F(0, 0.0, ρ)) ≈ sin(ρ) atol = 1.0e-8
            @test abs(imag(F(0, 0.0, ρ))) ≤ 1.0e-8
        end
    end

    @testset "D⁺ and D⁻ symmetry" begin
        for (ℓ, η) in ((0, 0.5), (1, 1.2), (2, -0.7))
            @test conj(D⁺(ℓ, η)) ≈ D⁻(ℓ, η) atol = 1.0e-10
        end
    end

    @testset "H⁺ and H⁻ symmetry" begin
        for (ℓ, η, ρ) in ((0, 0.5, 1.0), (1, 1.0, 2.0), (1.5, -0.3, 3.5))
            @test conj(H⁺(ℓ, η, ρ)) ≈ H⁻(ℓ, η, ρ) atol = 1.0e-10
        end
    end

    @testset "F_imag and G combination" begin
        for (ℓ, η, ρ) in ((0, 0.5, 1.0), (1, 1.0, 2.0), (1.5, -0.7, 3.5))
            @test F_imag(ℓ, η, ρ) ≈ F(ℓ, η, ρ) atol = 1.0e-10
            @test G(ℓ, η, ρ) ≈ real(G(ℓ, η, ρ)) atol = 1.0e-10  # G should be real
        end
    end

    @testset "M_regularized" begin
        α, β, γ = 2.0, 3.0, 4.0
        @test M_regularized(α, β, γ) ≈ HypergeometricFunctions._₁F₁(α, β, γ) / gamma(β) atol = 1.0e-10
    end

    @testset "Φ modified function" begin
        # consistency check vs definition: imaginary part small for real inputs
        for (ℓ, η, ρ) in ((0, 0.5, 1.0), (1, 1.0, 2.0))
            val = Φ(ℓ, η, ρ)
            @test abs(imag(val)) ≤ 1.0e-10
        end
    end

    @testset "Gamow factor w" begin
        # integer ℓ
        @test w(2, 1.0) == (1 + 0^2 / 1^2) * (1 + 1^2 / 1^2) * (1 + 2^2 / 1^2)
        # half-integer ℓ=1/2
        @test w(0.5, 2.0) ≈ (1 + 0.5^2 / 2.0^2) atol = 1.0e-10
    end


    @testset "Real‐valued F/G overloads" begin
        for (ℓ, η, ρ) in ((0.0, 0.5, 1.2), (1.0, -1.3, 2.5), (2.5, 0.7, 3.7))
            f_scalar = real(F(complex(ℓ), complex(η), complex(ρ)))
            @test isa(F(ℓ, η, ρ), Float64)
            @test F(ℓ, η, ρ) ≈ f_scalar atol = 1.0e-10

            g_scalar = real(G(complex(ℓ), complex(η), complex(ρ)))
            @test isa(G(ℓ, η, ρ), Float64)
            @test G(ℓ, η, ρ) ≈ g_scalar atol = 1.0e-10
        end
    end

    @testset "w_plus / w_minus symmetry" begin
        for (ℓ, η) in ((0, 0.7), (1.5, -1.2), (2.0, 0.3))
            @test conj(w_plus(ℓ, η)) ≈ w_minus(ℓ, η) atol = 1.0e-10
        end
    end

    @testset "h_plus / h_minus symmetry" begin
        for (ℓ, η) in ((0, 0.7), (1.5, 1.2), (2.0, -0.4))
            @test conj(h_plus(ℓ, η)) ≈ h_minus(ℓ, η) atol = 1.0e-10
        end
    end

    @testset "g real‐valued" begin
        for (ℓ, η) in ((0, 0.7), (1.5, 1.2), (2.0, -0.4))
            val = g(ℓ, η)
            @test isa(val, Real)
            @test abs(imag(val)) ≤ 1.0e-10
        end
    end


    ℓs = [0.0, 1.0, 2.5]
    ηs = [0.1, 1.0, 5.0]
    ρs = [0.5, 2.0, 10.0]

    @testset "Φ_dot numerical derivative" begin
        for ℓ in ℓs, η in ηs, ρ in ρs
            h = 1.0e-6
            # Central difference should be close to forward/backward difference for smooth functions
            central = Φ_dot(ℓ, η, ρ; h = h)
            forward = (Φ(ℓ + h, η, ρ) - Φ(ℓ, η, ρ)) / h
            backward = (Φ(ℓ, η, ρ) - Φ(ℓ - h, η, ρ)) / h
            @test isfinite(central)
            @test abs(central - 0.5 * (forward + backward)) < 1.0e-4 * max(1, abs(central))
        end
    end

    @testset "F_dot numerical derivative" begin
        for ℓ in ℓs, η in ηs, ρ in ρs
            h = 1.0e-6
            central = F_dot(ℓ, η, ρ; h = h)
            forward = (F(ℓ + h, η, ρ) - F(ℓ, η, ρ)) / h
            backward = (F(ℓ, η, ρ) - F(ℓ - h, η, ρ)) / h
            @test isfinite(central)
            @test abs(central - 0.5 * (forward + backward)) < 1.0e-4 * max(1, abs(central))
        end
    end


    @testset "Φ_dot and F_dot handle keyword argument h" begin
        ℓ, η, ρ = 1.0, 1.0, 1.0
        val1 = Φ_dot(ℓ, η, ρ)
        val2 = Φ_dot(ℓ, η, ρ; h = 1.0e-5)
        @test isfinite(val1)
        @test isfinite(val2)
        val3 = F_dot(ℓ, η, ρ)
        val4 = F_dot(ℓ, η, ρ; h = 1.0e-5)
        @test isfinite(val3)
        @test isfinite(val4)
    end
end
