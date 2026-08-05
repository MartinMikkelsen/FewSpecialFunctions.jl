using Test, DelimitedFiles, FewSpecialFunctions, SpecialFunctions

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
