using Test
using FewSpecialFunctions
using ForwardDiff
using SpecialFunctions
import FewSpecialFunctions: dawson

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
            ForwardDiff.derivative(x -> debye_function(2.0, 2.0, x), 0.5),
            fdiff(x -> debye_function(2.0, 2.0, x), 0.5);
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

    @testset "Dawson integral" begin
        for x in (-2.0, -0.1, 0.1, 1.0, 10.0)
            expected = muladd(-2x, dawson(x), 1.0)
            @test isapprox(ForwardDiff.derivative(dawson, x), expected; rtol = 1.0e-12)
        end
        @test dawson(ForwardDiff.Dual{Nothing}(1.0, 1.0)) isa ForwardDiff.Dual
    end

    @testset "Voigt function" begin
        x, y = 1.0, 0.5
        @test isapprox(ForwardDiff.derivative(t -> voigt(t, y), x), fdiff(t -> voigt(t, y), x); rtol = 1.0e-5)
        @test isapprox(ForwardDiff.derivative(t -> voigt(x, t), y), fdiff(t -> voigt(x, t), y); rtol = 1.0e-5)

        xdual = ForwardDiff.Dual{Nothing}(x, 1.0, 0.0)
        ydual = ForwardDiff.Dual{Nothing}(y, 0.0, 1.0)
        zdual = voigt(xdual, ydual)
        @test zdual isa ForwardDiff.Dual
        @test isapprox(ForwardDiff.value(zdual), voigt(x, y); rtol = 1.0e-12)
        @test isapprox(ForwardDiff.partials(zdual)[1], fdiff(t -> voigt(t, y), x); rtol = 1.0e-5)
        @test isapprox(ForwardDiff.partials(zdual)[2], fdiff(t -> voigt(x, t), y); rtol = 1.0e-5)

        xunit = ForwardDiff.Dual{Nothing}(1.0, 1.0)
        yunit = ForwardDiff.Dual{Nothing}(0.5, 1.0)
        zunit = voigt(xunit, yunit)
        @test zunit isa ForwardDiff.Dual
        @test isapprox(ForwardDiff.value(zunit), voigt(1.0, 0.5); rtol = 1.0e-12)
        @test isapprox(ForwardDiff.partials(zunit)[1], fdiff(t -> voigt(t, 0.5), 1.0) + fdiff(t -> voigt(1.0, t), 0.5); rtol = 1.0e-5)
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

    @testset "Bose–Einstein family" begin
        for k in (-4.5, -3.0, -0.5, 0.0, 0.5, 1.5, 5.5, 19.5, 20.0), x in (-2.5, -0.1)
            # Independent derivative of the fugacity expansion, evaluated at
            # high precision; 3000 terms leave a negligible exponential tail.
            reference = Float64(sum(exp(big(n) * big(x)) / big(n)^big(k) for n in 1:3000))
            @test ForwardDiff.derivative(t -> BoseEinsteinIntegralNorm(k, t), x) ≈ reference rtol = 5.0e-15
            if k > -1
                @test ForwardDiff.derivative(t -> BoseEinsteinIntegral(k, t), x) ≈ gamma(k + 1) * reference rtol = 5.0e-15
            end
        end
        for k in (-4.5, 0.5, 1.5)
            reference = Float64(sum(exp(-big(n)) / big(n)^(big(k) - 1) for n in 1:300))
            second = ForwardDiff.derivative(x -> ForwardDiff.derivative(t -> BoseEinsteinIntegralNorm(k, t), x), -1.0)
            @test second ≈ reference rtol = 5.0e-15
            if k > -1
                second_scaled = ForwardDiff.derivative(x -> ForwardDiff.derivative(t -> BoseEinsteinIntegral(k, t), x), -1.0)
                @test second_scaled ≈ gamma(k + 1) * reference rtol = 5.0e-15
            end
        end
        @test ForwardDiff.derivative(x -> BoseEinsteinIntegralNorm(1.5f0, x), -1.0f0) isa Float32
        for (k, x) in ((35.0f0, -10.0f0), (40.0f0, -30.0f0))
            reference = Float32(gamma(BigFloat(k) + 1) * BoseEinsteinIntegralNorm(BigFloat(k) - 1, BigFloat(x)))
            @test ForwardDiff.derivative(t -> BoseEinsteinIntegral(k, t), x) ≈ reference rtol = 2eps(Float32)
        end
        @test ForwardDiff.derivative(x -> BoseEinsteinIntegralNorm(1.5, x), big"-1") isa BigFloat
        @test ForwardDiff.derivative(x -> BoseEinsteinIntegralNorm(1.5, x), 0.0) ≈ zeta(1.5)
        @test ForwardDiff.derivative(x -> BoseEinsteinIntegralNorm(0.5, x), 0.0) == Inf
        @test ForwardDiff.derivative(x -> BoseEinsteinIntegral(200, x), -1000.0) ≈ Float64(gamma(big(201)) * exp(big(-1000))) rtol = 3.0e-13
        for f in (BoseEinsteinIntegral, BoseEinsteinIntegralNorm)
            @test ForwardDiff.derivative(x -> f(UInt(0), x), -0.5) ≈ inv(expm1(0.5))
            @test ForwardDiff.derivative(x -> f(UInt(0), x), -0.5f0) isa Float32
            @test ForwardDiff.derivative(x -> ForwardDiff.derivative(t -> f(UInt(1), t), x), -0.5) ≈ inv(expm1(0.5))
            @test_throws DomainError ForwardDiff.derivative(k -> f(k, -1.0), 0.5)
            @test_throws DomainError ForwardDiff.gradient(v -> f(v[1], v[2]), [0.5, -1.0])
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
        # Independent Gaussian and Bessel identities at the reported failures.
        for x in (-10.0, 10.0)
            @test ForwardDiff.derivative(t -> U(-0.5, t), x) ≈ -x / 2 * exp(-x^2 / 4) rtol = 2.0e-13
        end
        w10 = sqrt(10pi) / 2^(5 / 4) * (besselj(-0.25, 25) - besselj(0.25, 25))
        dw10 = -10sqrt(10pi) / 2^(9 / 4) * (besselj(-0.75, 25) + besselj(0.75, 25))
        @test ForwardDiff.derivative(x -> W(0.0, x), 10.0) ≈ dw10 rtol = 2.0e-12
        @test ForwardDiff.derivative(x -> dW(0.0, x), 10.0) ≈ -25w10 rtol = 2.0e-12

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

        # Clausen n=1 uses finite-difference derivative (no closed form)
        @testset "Clausen n=1 ForwardDiff" begin
            θ0 = π / 3
            fd = fdiff(x -> Clausen(1, x), θ0)
            ad = ForwardDiff.derivative(x -> Clausen(1, x), θ0)
            @test isapprox(ad, fd; rtol = 1.0e-5)
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
