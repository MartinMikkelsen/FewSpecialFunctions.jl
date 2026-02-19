using Test
using FewSpecialFunctions
using ForwardDiff
using SpecialFunctions

fdiff(f, x; h = 1.0e-6) = (f(x + h) - f(x - h)) / (2h)

@testset "ForwardDiff extension" begin
    @testset "Coulomb family" begin
        @test isapprox(ForwardDiff.derivative(x -> η(x), 2.0), -1 / (2 * 2.0^(3 / 2)); atol = 1.0e-8)
        @test isapprox(ForwardDiff.derivative(x -> η(3.0, x), 1.7), -1 / (3.0 * 1.7^2); atol = 1.0e-8)
        @test isapprox(ForwardDiff.derivative(x -> η(x, 2.0), 3.0), -1 / (3.0^2 * 2.0); atol = 1.0e-8)

        @test isfinite(ForwardDiff.derivative(x -> real(C(0.5, x)), 0.2))
        @test isfinite(ForwardDiff.derivative(x -> real(D⁺(0.5, x)), 0.2))
        @test isfinite(ForwardDiff.derivative(x -> real(D⁻(0.5, x)), 0.2))

        @test isapprox(ForwardDiff.derivative(x -> θ(0.5, 0.2, x), 1.4), 1 - 0.2 / 1.4; atol = 1.0e-7)
        @test isapprox(ForwardDiff.derivative(x -> θ(1.0, 0.5, x), 2.0), 1 - 0.5 / 2.0; atol = 1.0e-7)
        @test isapprox(ForwardDiff.derivative(x -> θ(2.0, 0.1, x), 3.0), 1 - 0.1 / 3.0; atol = 1.0e-7)

        @test isapprox(ForwardDiff.derivative(x -> F(0.0, 0.2, x), 1.0), fdiff(x -> F(0.0, 0.2, x), 1.0); atol = 1.0e-5, rtol = 1.0e-4)
        @test isapprox(ForwardDiff.derivative(x -> F(1.0, 0.3, x), 2.0), fdiff(x -> F(1.0, 0.3, x), 2.0); atol = 1.0e-5, rtol = 1.0e-4)
        @test isapprox(ForwardDiff.derivative(x -> G(0.0, 0.2, x), 1.5), fdiff(x -> G(0.0, 0.2, x), 1.5); atol = 1.0e-5, rtol = 1.0e-4)
        @test isapprox(ForwardDiff.derivative(x -> F(0.0, x, 1.5), 0.2), fdiff(x -> F(0.0, x, 1.5), 0.2); atol = 1.0e-5, rtol = 1.0e-4)

        @test isfinite(ForwardDiff.derivative(x -> real(H⁺(0.5, 0.2, x)), 1.0))
        @test isfinite(ForwardDiff.derivative(x -> real(H⁻(0.5, 0.2, x)), 1.0))
        @test isfinite(ForwardDiff.derivative(x -> real(F_imag(0.5, 0.2, x)), 1.0))
        @test isfinite(ForwardDiff.derivative(x -> M_regularized(0.5, 1.2, x), 0.4))
        @test isfinite(ForwardDiff.derivative(x -> real(Φ(0.5, 0.2, x)), 1.0))
        @test isfinite(ForwardDiff.derivative(x -> w(2, x), 1.5))
    end

    @testset "Debye family" begin
        d1 = ForwardDiff.derivative(x -> debye_function(2.0, 1.0, x), 1.2)
        r1 = fdiff(x -> debye_function(2.0, 1.0, x), 1.2)
        @test isapprox(d1, r1; atol = 1.0e-5, rtol = 1.0e-4)

        d2 = ForwardDiff.derivative(x -> debye_function(x, 1.0), 1.2)
        r2 = fdiff(x -> debye_function(x, 1.0), 1.2)
        @test isapprox(d2, r2; atol = 1.0e-5, rtol = 1.0e-4)

        # Different orders and β values
        @test isapprox(
            ForwardDiff.derivative(x -> debye_function(3.0, 1.0, x), 2.0),
            fdiff(x -> debye_function(3.0, 1.0, x), 2.0);
            atol = 1.0e-5,
            rtol = 1.0e-4,
        )
        @test isapprox(
            ForwardDiff.derivative(x -> debye_function(1.0, 2.0, x), 0.5),
            fdiff(x -> debye_function(1.0, 2.0, x), 0.5);
            atol = 1.0e-5,
            rtol = 1.0e-4,
        )
        @test isapprox(
            ForwardDiff.derivative(x -> debye_function(4.0, 1.0, x), 3.0),
            fdiff(x -> debye_function(4.0, 1.0, x), 3.0);
            atol = 1.0e-5,
            rtol = 1.0e-4,
        )
    end

    @testset "Fresnel family" begin
        x0 = 0.3
        @test isapprox(ForwardDiff.derivative(FresnelC, x0), cos((π / 2) * x0^2); atol = 1.0e-10)
        @test isapprox(ForwardDiff.derivative(FresnelS, x0), sin((π / 2) * x0^2); atol = 1.0e-10)
        @test isapprox(ForwardDiff.derivative(x -> real(FresnelE(x)), x0), real(exp(im * (π / 2) * x0^2)); atol = 1.0e-10)
        @test isapprox(ForwardDiff.derivative(x -> imag(FresnelE(x)), x0), imag(exp(im * (π / 2) * x0^2)); atol = 1.0e-10)

        # Analytic derivatives hold at multiple evaluation points
        for x1 in [0.1, 0.5, 1.0, 1.5, 2.0]
            @test isapprox(ForwardDiff.derivative(FresnelC, x1), cos((π / 2) * x1^2); atol = 1.0e-10)
            @test isapprox(ForwardDiff.derivative(FresnelS, x1), sin((π / 2) * x1^2); atol = 1.0e-10)
        end

        Cx, Sx, Ex = fresnel(ForwardDiff.Dual{Nothing}(x0, 1.0))
        @test Cx isa ForwardDiff.Dual
        @test Sx isa ForwardDiff.Dual
        @test Ex isa Complex
        @test real(Ex) isa ForwardDiff.Dual
        @test imag(Ex) isa ForwardDiff.Dual
    end

    @testset "Clausen family" begin
        d = ForwardDiff.derivative(x -> Clausen(2, x), 1.1)
        r = fdiff(x -> Clausen(2, x), 1.1)
        @test isapprox(d, r; atol = 1.0e-5, rtol = 1.0e-4)

        # Recursive analytic identity: dCl_n/dθ = Cl_{n-1}(θ) (n even), -Cl_{n-1}(θ) (n odd ≥ 3)
        for θ0 in [0.5, 1.1, 2.3]
            @test isapprox(ForwardDiff.derivative(x -> Clausen(2, x), θ0), Clausen(1, θ0); atol = 1.0e-5, rtol = 1.0e-4)
            @test isapprox(ForwardDiff.derivative(x -> Clausen(3, x), θ0), -Clausen(2, θ0); atol = 1.0e-5, rtol = 1.0e-4)
            @test isapprox(ForwardDiff.derivative(x -> Clausen(4, x), θ0), Clausen(3, θ0); atol = 1.0e-5, rtol = 1.0e-4)
        end
    end

    @testset "Fermi-Dirac family" begin
        d = ForwardDiff.derivative(x -> FermiDiracIntegral(1.5, x), 0.2)
        r = fdiff(x -> FermiDiracIntegral(1.5, x), 0.2)
        @test isapprox(d, r; atol = 1.0e-5, rtol = 1.0e-4)
        @test isfinite(ForwardDiff.derivative(x -> FermiDiracIntegralNorm(1.5, x), 0.2))

        # Multiple half-integer orders and evaluation points
        for j in [-0.5, 0.5, 2.5]
            for x0 in [-1.0, 0.2, 1.5]
                @test isapprox(
                    ForwardDiff.derivative(x -> FermiDiracIntegral(j, x), x0),
                    fdiff(x -> FermiDiracIntegral(j, x), x0);
                    atol = 1.0e-5,
                    rtol = 1.0e-4,
                )
            end
        end
        @test isapprox(
            ForwardDiff.derivative(x -> FermiDiracIntegralNorm(1.5, x), 0.2),
            fdiff(x -> FermiDiracIntegralNorm(1.5, x), 0.2);
            atol = 1.0e-5,
            rtol = 1.0e-4,
        )
    end

    @testset "MarcumQ family" begin
        b0 = 1.2
        @test isapprox(ForwardDiff.derivative(x -> MarcumQ(1.0, 2.0, x), b0), dQdb(1.0, 2.0, b0); atol = 1.0e-5, rtol = 1.0e-4)
        @test isfinite(ForwardDiff.derivative(x -> dQdb(1.0, 2.0, x), b0))

        # Test d/db at several (M, a, b) combinations
        for (M, a, b) in [(2.0, 1.5, 3.0), (3.0, 0.5, 2.0), (1.0, 3.0, 4.0)]
            @test isapprox(ForwardDiff.derivative(x -> MarcumQ(M, a, x), b), dQdb(M, a, b); atol = 1.0e-5, rtol = 1.0e-4)
        end

        # Derivative w.r.t. a compared with finite difference
        @test isapprox(
            ForwardDiff.derivative(a -> MarcumQ(1.0, a, 1.2), 2.0),
            fdiff(a -> MarcumQ(1.0, a, 1.2), 2.0);
            atol = 1.0e-5,
            rtol = 1.0e-4,
        )

        # 2-arg convenience form
        @test isapprox(ForwardDiff.derivative(x -> MarcumQ(1.0, x), 1.5), dQdb(1.0, 1.5); atol = 1.0e-5, rtol = 1.0e-4)
    end

    @testset "Parabolic cylinder family" begin
        x0 = 0.7
        @test isapprox(ForwardDiff.derivative(x -> U(0.2, x), x0), dU(0.2, x0); atol = 1.0e-5, rtol = 1.0e-4)
        @test isapprox(ForwardDiff.derivative(x -> V(0.2, x), x0), dV(0.2, x0); atol = 1.0e-5, rtol = 1.0e-4)
        @test isapprox(ForwardDiff.derivative(x -> W(0.2, x), x0), dW(0.2, x0); atol = 1.0e-5, rtol = 1.0e-4)

        @test isfinite(ForwardDiff.derivative(x -> dU(0.2, x), x0))
        @test isfinite(ForwardDiff.derivative(x -> dV(0.2, x), x0))
        @test isfinite(ForwardDiff.derivative(x -> dW(0.2, x), x0))

        # Test at negative x and different a values
        for (a, x1) in [(0.2, -0.7), (1.5, 1.0), (0.5, 1.2)]
            @test isapprox(ForwardDiff.derivative(x -> U(a, x), x1), dU(a, x1); atol = 1.0e-5, rtol = 1.0e-4)
            @test isapprox(ForwardDiff.derivative(x -> V(a, x), x1), dV(a, x1); atol = 1.0e-5, rtol = 1.0e-4)
        end

        # Second derivative from the parabolic cylinder ODE: d(dU)/dx = (x²/4 + a) * U(a, x)
        for (a, x1) in [(0.2, 0.7), (0.5, 1.2), (1.5, -0.8)]
            d2U = ForwardDiff.derivative(x -> dU(a, x), x1)
            @test isapprox(d2U, (x1^2 / 4 + a) * U(a, x1); atol = 1.0e-5, rtol = 1.0e-4)
            d2V = ForwardDiff.derivative(x -> dV(a, x), x1)
            @test isapprox(d2V, (x1^2 / 4 + a) * V(a, x1); atol = 1.0e-5, rtol = 1.0e-4)
        end
    end
end
