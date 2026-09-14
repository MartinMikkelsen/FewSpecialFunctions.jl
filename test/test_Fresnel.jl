using Test, DelimitedFiles, FewSpecialFunctions, SpecialFunctions
using QuadGK: quadgk

@testset "Fresnel" begin

    @testset "Fresnel values" begin
        # read 5 columns: x, fsin, fcos, fexp_real, fexp_imag
        data = open(readdlm, joinpath(@__DIR__, "data", "FresnelF.txt"))

        for r in 1:size(data, 1)
            x = data[r, 1]
            fsin_ref = data[r, 2]
            fcos_ref = data[r, 3]
            # Ensure complex parsing is robust
            fexp_str = string(data[r, 4])
            fexp_ref = parse(ComplexF64, replace(fexp_str, "i" => "im"))

            @test isapprox(FewSpecialFunctions.FresnelS(x), fsin_ref; rtol = 1.0e-2)
            @test isapprox(FewSpecialFunctions.FresnelC(x), fcos_ref; rtol = 1.0e-2)
        end
    end

end
@testset "fresnel function" begin
    # Test fresnel function directly
    x_values = [0.0, 1.0, 2.0, -1.5]

    for x in x_values
        C, S, E = FewSpecialFunctions.fresnel(x)

        # Test individual components match wrapper functions
        @test C ≈ FewSpecialFunctions.FresnelC(x)
        @test S ≈ FewSpecialFunctions.FresnelS(x)
        @test E ≈ FewSpecialFunctions.FresnelE(x)

        # Test that E = C + i*S relationship holds
        @test E ≈ C + im * S
    end

    # Test special case z=0
    C, S, E = FewSpecialFunctions.fresnel(0.0)
    @test C ≈ 0.0
    @test S ≈ 0.0
    @test E ≈ 0.0

    # Test complex input
    z = 1.0 + 1.0im
    C, S, E = FewSpecialFunctions.fresnel(z)
    @test E ≈ C + im * S

    # Test against known values
    # At z = 1, C(1) ≈ 0.7798934, S(1) ≈ 0.4382591
    C, S, E = FewSpecialFunctions.fresnel(1.0)
    @test isapprox(C, 0.7798934, rtol = 1.0e-6)
    @test isapprox(S, 0.4382591, rtol = 1.0e-6)
end

@testset "Fresnel edge and special cases" begin

    # Test at purely imaginary input
    z = 2.0im
    C, S, E = FewSpecialFunctions.fresnel(z)
    # E should equal C + im*S
    @test E ≈ C + im * S

    # Test at a small value (Taylor expansion regime)
    z = 1.0e-8
    C, S, E = FewSpecialFunctions.fresnel(z)
    @test isapprox(C, z, atol = 1.0e-8)
    @test isapprox(S, (π / 6) * z^3, atol = 1.0e-24)

    # Test at a negative real value
    z = -2.0
    C, S, E = FewSpecialFunctions.fresnel(z)
    # Fresnel integrals are odd/even functions:
    @test isapprox(C, -FewSpecialFunctions.FresnelC(2.0), rtol = 1.0e-6)
    @test isapprox(S, -FewSpecialFunctions.FresnelS(2.0), rtol = 1.0e-6)

end

@testset "Fresnel complex values" begin
    z = 1.0 + 1.0im
    C, S, E = fresnel(z)

    @test C isa ComplexF64
    @test S isa ComplexF64
    @test isapprox(C, 2.5557937781024376 + 2.5557937781024376im; rtol = 1.0e-14)
    @test isapprox(S, -2.0618882191948393 + 2.0618882191948393im; rtol = 1.0e-14)
    @test E ≈ C + im * S
    @test FresnelC(z) ≈ C
    @test FresnelS(z) ≈ S
    @test FresnelE(z) ≈ E
end

@testset "Fresnel complex intermediate value" begin
    z = 2.0 + 3.0im
    C, S, E = fresnel(z)

    @test isapprox(C, -3.788100200182899e6 + 5.815899102940467e6im; rtol = 1.0e-13)
    @test isapprox(S, -5.815898602940467e6 - 3.788100700182899e6im; rtol = 1.0e-13)
    @test E ≈ C + im * S
end

@testset "Fresnel complex cancellation and sectors" begin
    # Independent references from the defining integrals, evaluated with
    # 320-bit BigFloat quadrature. E must not be reconstructed from C and S
    # where the components cancel many digits.
    @test FresnelE(4 + 4im) ≈ 0.5 + 0.5im rtol = 1.0e-14
    @test FresnelE(8 + 8im) ≈ 0.5 + 0.5im rtol = 1.0e-14
    Cref = 2.521551371184569e13 - 3.107276588837978e14im
    Sref = 3.107276588837983e14 + 2.521551371184519e13im
    @test FresnelC(1 + 12im) ≈ Cref rtol = 1.0e-13
    @test FresnelS(1 + 12im) ≈ Sref rtol = 1.0e-13
    @test FresnelC(1 - 12im) ≈ conj(Cref) rtol = 1.0e-13
    @test FresnelS(1 - 12im) ≈ conj(Sref) rtol = 1.0e-13

    for z in (1.0 + 12im, 12.0 + im, 2.0 + 3im, 4.0 + 4im)
        C, S, E = fresnel(z)
        @test FresnelC(conj(z)) ≈ conj(C) rtol = 1.0e-13
        @test FresnelS(conj(z)) ≈ conj(S) rtol = 1.0e-13
        @test FresnelC(im * z) ≈ im * C rtol = 1.0e-13
        @test FresnelS(im * z) ≈ -im * S rtol = 1.0e-13
        @test FresnelE(-z) ≈ -E rtol = 1.0e-13
        # These identities remain well conditioned when E is small.
        Eminus = conj(FresnelE(conj(z)))
        @test E + Eminus ≈ 2C rtol = 1.0e-13
        @test E - Eminus ≈ 2im * S rtol = 1.0e-13
    end

    for T in (Float32, Float64, BigFloat)
        vals = fresnel(complex(T(1), T(1)))
        @test all(v -> v isa Complex{T}, vals)
        @test vals[3] ≈ complex(T(0.49390555890759856), T(0.49390555890759856)) rtol = 20eps(T) + eps(Float64)
    end

    # The old complex series/Miller and Miller/asymptotic switches.
    for radius in sqrt.((6.9, 25 + 2precision(Float64)))
        z = radius * cis(0.7)
        for f in (FresnelC, FresnelS, FresnelE)
            @test f((1 - 1.0e-10) * z) ≈ f((1 + 1.0e-10) * z) rtol = 1.0e-7
            @test f((1 - 1.0e-10) * conj(z)) ≈ f((1 + 1.0e-10) * conj(z)) rtol = 1.0e-7
        end
    end
end

@testset "BigFloat Fresnel conjugation across numerical regions" begin
    setprecision(BigFloat, 128) do
        # Series, Miller recurrence, and the imaginary-dominant asymptotic
        # sector. Integrate along a straight complex segment independently.
        for (x, y) in ((1, 1), (2, 3), (1 // 10, 18))
            z = complex(BigFloat(x), BigFloat(y))
            Cref, _ = quadgk(t -> z * cos(BigFloat(π) * (z * t)^2 / 2), zero(BigFloat), one(BigFloat); rtol = big"1e-32")
            Sref, _ = quadgk(t -> z * sin(BigFloat(π) * (z * t)^2 / 2), zero(BigFloat), one(BigFloat); rtol = big"1e-32")
            C, S, E = fresnel(conj(z))
            @test C ≈ conj(Cref) rtol = big"1e-30"
            @test S ≈ conj(Sref) rtol = big"1e-30"
            @test E ≈ conj(Cref) + im * conj(Sref) rtol = big"1e-30"
        end
    end
end

@testset "Fresnel dispatch" begin
    Cfloat, Sfloat, Efloat = fresnel(1.0)
    @test Cfloat isa Float64
    @test Sfloat isa Float64
    @test Efloat isa ComplexF64
    @test isapprox(Sfloat, 0.4382591473903548; rtol = 1.0e-15)
    @test FresnelS(1.0) == Sfloat

    Cint, Sint, Eint = fresnel(1)
    @test Cint isa Float64
    @test Sint isa Float64
    @test Eint isa ComplexF64
    @test (Cint, Sint, Eint) == (Cfloat, Sfloat, Efloat)

    zfloat = 1.0 + 1.0im
    Ccomplex, Scomplex, Ecomplex = fresnel(zfloat)
    @test Ccomplex isa ComplexF64
    @test Scomplex isa ComplexF64
    @test Ecomplex isa ComplexF64
    @test FresnelS(zfloat) == Scomplex

    Ccomplexint, Scomplexint, Ecomplexint = fresnel(1 + im)
    @test Ccomplexint isa ComplexF64
    @test Scomplexint isa ComplexF64
    @test Ecomplexint isa ComplexF64
    @test (Ccomplexint, Scomplexint, Ecomplexint) == (Ccomplex, Scomplex, Ecomplex)

    Creal, Sreal, Ereal = fresnel(1 // 2)
    @test Creal isa Float64
    @test Sreal isa Float64
    @test Ereal isa ComplexF64
    @test (Creal, Sreal, Ereal) == fresnel(0.5)

    zrational = complex(1 // 2, 1 // 3)
    Cgeneric, Sgeneric, Egeneric = fresnel(zrational)
    @test Cgeneric isa ComplexF64
    @test Sgeneric isa ComplexF64
    @test Egeneric isa ComplexF64
    @test (Cgeneric, Sgeneric, Egeneric) == fresnel(complex(0.5, 1 / 3))
end

@testset "Fresnel asymptotic and complex branches" begin
    # Real intermediate region: the error-function evaluation path.
    Cintermediate, Sintermediate, Eintermediate = fresnel(3.0)
    wintermediate = (sqrt(π) / 2) * (1 - im) * 3.0
    Eintermediate_ref = (1 + im) / 2 * erf(wintermediate)
    @test Eintermediate ≈ Eintermediate_ref
    @test Cintermediate + im * Sintermediate == Eintermediate

    # Real large-argument region: the asymptotic evaluation path.
    Casymptotic, Sasymptotic, Easymptotic = fresnel(10.0)
    @test isapprox(Casymptotic, 0.49989869420551575; rtol = 1.0e-14)
    @test isapprox(Sasymptotic, 0.46816997858488224; rtol = 1.0e-14)
    @test Casymptotic + im * Sasymptotic == Easymptotic

    # Complex large-argument region with |imag(z)| > |real(z)|.
    zasymptotic = 1.0 + 12.0im
    Casymptotic_complex, Sasymptotic_complex, Easymptotic_complex = fresnel(zasymptotic)
    wasymptotic = (sqrt(π) / 2) * (1 - im) * zasymptotic
    Easymptotic_ref = (1 + im) / 2 * erf(wasymptotic)
    @test isapprox(Easymptotic_complex, Easymptotic_ref; rtol = 1.0e-13)
    @test Casymptotic_complex + im * Sasymptotic_complex == Easymptotic_complex

    # Complex dispatch branches for zero imaginary and negative real parts.
    Creal_complex, Sreal_complex, Ereal_complex = fresnel(3.0 + 0.0im)
    @test (Creal_complex, Sreal_complex, Ereal_complex) ==
        (complex(Cintermediate), complex(Sintermediate), Eintermediate)

    Cnegative, Snegative, Enegative = fresnel(-1.0 + 1.0im)
    Cpositive, Spositive, Epositive = fresnel(1.0 - 1.0im)
    @test Cnegative == -Cpositive
    @test Snegative == -Spositive
    @test Enegative == -Epositive
end

@testset "Fresnel helper functions" begin
    # The real error-function decomposition reproduces known C(3) and S(3).
    Cerf, Serf = FewSpecialFunctions._fresnel_erf_real(3.0)
    @test isapprox(Cerf, 0.6057207892976857; rtol = 1.0e-14)
    @test isapprox(Serf, 0.496312998967375; rtol = 1.0e-14)

    # The asymptotic threshold is precision-aware and never below five.
    @test isapprox(
        FewSpecialFunctions._fresnel_asymptotic_start(Float64),
        5.790209015886026;
        rtol = 1.0e-15,
    )
    @test FewSpecialFunctions._fresnel_asymptotic_start(BigFloat) > BigFloat(5)
end
