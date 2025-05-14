using DelimitedFiles
using SpecialFunctions

# some tests from https://github.com/Expander/ClausenFunctions.jl

function complex_approx(a::ComplexF64, b::ComplexF64; atol=1e-10)
    return isapprox(real(a), real(b), atol=atol) && isapprox(imag(a), imag(b), atol=atol)
end

@testset "Clausen1" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "Cl1.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(1,x;m=500000) ≈ expected rtol=1e-11
        @test FewSpecialFunctions.Clausen(1,-x) ≈ expected rtol=1e-12
        @test FewSpecialFunctions.Clausen(1,x-2pi) ≈ expected rtol=1e-12
        @test FewSpecialFunctions.Clausen(1,x+2pi) ≈ expected rtol=1e-12

    end

    @test FewSpecialFunctions.Clausen(1,0.5) ≈ 0.70358563513784466 rtol=1e-14
    @test FewSpecialFunctions.Clausen(1,0.5;N=20) ≈ 0.70358563513784466 rtol=1e-14
    @test FewSpecialFunctions.Clausen(1,1.0) ≈ 0.042019505825368962 rtol=1e-14
    @test FewSpecialFunctions.Clausen(1,2.0) ≈ -0.52054343429085363 rtol=1e-14

end

@testset "Clausen2" begin

    data2 = open(readdlm, joinpath(@__DIR__, "data", "Cl2.txt"))

    for r in 1:size(data2, 1)
        row      = data2[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(2,x;N=20,m=500000) ≈ expected rtol=1e-9

    end

end

@testset "Clausen3" begin

    data3 = open(readdlm, joinpath(@__DIR__, "data", "Cl3.txt"))

    for r in 1:size(data3, 1)
        row      = data3[r, :]
        x        = row[1]
        expected = row[2]
        
        @test FewSpecialFunctions.Clausen(3,x;N=10,m=500000) ≈ expected rtol=1e-9
        @test FewSpecialFunctions.Clausen(3,x;N=20,m=500000) ≈ expected rtol=1e-9

    end

end

@testset "Clausen4" begin

    data4 = open(readdlm, joinpath(@__DIR__, "data", "Cl4.txt"))

    for r in 1:size(data4, 1)
        row      = data4[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(4,x;N=10,m=500000) ≈ expected rtol=1e-8
        @test FewSpecialFunctions.Clausen(4,x;N=20,m=500000) ≈ expected rtol=1e-8
    end
end

@testset "Clausen5" begin

    data5 = open(readdlm, joinpath(@__DIR__, "data", "Cl5.txt"))

    for r in 1:size(data5, 1)
        row      = data5[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(5,x;N=10,m=500000) ≈ expected rtol=1e-8
        @test FewSpecialFunctions.Clausen(5,x;N=20,m=500000) ≈ expected rtol=1e-8
    end
end

@testset "Clausen6" begin

    data6 = open(readdlm, joinpath(@__DIR__, "data", "Cl6.txt"))

    for r in 1:size(data6, 1)
        row      = data6[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(6,x) ≈ expected rtol=1e-8
        @test FewSpecialFunctions.Clausen(6,x) ≈ expected rtol=1e-8
    end
end

@testset "f_n function" begin
    # Test with even n (should use sine)
    @test FewSpecialFunctions.f_n(2, 1, π/4) ≈ sin(π/4) rtol=1e-15
    @test FewSpecialFunctions.f_n(2, 2, π/3) ≈ sin(2π/3)/2^2 rtol=1e-15
    @test FewSpecialFunctions.f_n(4, 3, π/6) ≈ sin(3π/6)/3^4 rtol=1e-15
    @test FewSpecialFunctions.f_n(6, 5, 0.5) ≈ sin(5*0.5)/5^6 rtol=1e-15
    
    # Test with odd n (should use cosine)
    @test FewSpecialFunctions.f_n(1, 1, π/4) ≈ cos(π/4) rtol=1e-15
    @test FewSpecialFunctions.f_n(3, 2, π/3) ≈ cos(2π/3)/2^3 rtol=1e-15
    @test FewSpecialFunctions.f_n(5, 3, π/6) ≈ cos(3π/6)/3^5 rtol=1e-15
    @test FewSpecialFunctions.f_n(7, 5, 0.5) ≈ cos(5*0.5)/5^7 rtol=1e-15
    
    # Test special cases
    @test FewSpecialFunctions.f_n(2, 1, 0.0) == 0.0
    @test FewSpecialFunctions.f_n(1, 1, 0.0) == 1.0
    @test FewSpecialFunctions.f_n(2, 1, π) ≈ 0.0 atol=1e-15
    @test FewSpecialFunctions.f_n(1, 1, π) ≈ -1.0 rtol=1e-15
    
    # Test with larger k values
    @test FewSpecialFunctions.f_n(2, 10, 0.1) ≈ sin(10*0.1)/10^2 rtol=1e-14
    @test FewSpecialFunctions.f_n(3, 10, 0.1) ≈ cos(10*0.1)/10^3 rtol=1e-14
end

@testset "Testing F function dispatch" begin
    # Test inputs
    z = 1.0 + 2.0im
    θ = π / 4

    # Expected values for F1 to F6
    expected_F1 = Ci_complex(z * θ)
    expected_F2 = (θ * z * Ci_complex(z * θ) - sin(θ * z)) / z
    expected_F3 = -1 / (2 * z^2) * (θ^2 * z^2 * Ci_complex(z * θ) + cos(θ * z) - θ * z * sin(θ * z))
    expected_F4 = -1 / (6 * z^3) * (θ^3 * z^3 * Ci_complex(z * θ) + (2 - θ^2 * z^2) * sin(θ * z) + θ * z * cos(θ * z))
    expected_F5 = 1 / (24 * z^4) * (
        (θ^4 * z^4 * Ci_complex(z * θ)) +
        (θ * z * (2 - θ^2 * z^2) * sin(θ * z)) +
        ((θ^2 * z^2 - 6) * cos(θ * z))
    )
    expected_F6 = 1 / (120 * z^5) * (θ^5 * z^5 * Ci_complex(z * θ)+ θ*z*(θ^2*z^2-6)*cos(θ*z)-(θ^4*z^4-2*θ^2*z^2+24)*sin(θ*z))

    # Valid dispatch tests for F1 to F6
    @test complex_approx(F(1, z, θ), expected_F1)
    @test complex_approx(F(2, z, θ), expected_F2)
    @test complex_approx(F(3, z, θ), expected_F3)
    @test complex_approx(F(4, z, θ), expected_F4)
    @test complex_approx(F(5, z, θ), expected_F5)
    @test complex_approx(F(6, z, θ), expected_F6)


    # Test invalid n (less than 1 and greater than 6)
    @test_throws ArgumentError F(0, z, θ)
    @test_throws ArgumentError F(7, z, θ)

    # Test non-integer n (should throw MethodError)
    @test_throws MethodError F(2.5, z, θ)

    # Test negative n (should throw ArgumentError)
    @test_throws ArgumentError F(-1, z, θ)

end

@testset "Ci_complex Tests" begin
    # Test case 1: Negative real axis
    z2 = -1.0 + 0.0im
    iz2 = im * z2
    expected2 = -0.5 * (expint(iz2) + expint(-iz2)) + π * im
    @test complex_approx(Ci_complex(z2), expected2)


    # Test case 4: Positive real and imaginary parts
    z5 = 1.0 + 1.0im
    iz5 = im * z5
    expected5 = -0.5 * (expint(iz5) + expint(-iz5))
    @test complex_approx(Ci_complex(z5), expected5)

    # Test case 5: Negative real and positive imaginary parts
    z6 = -1.0 + 1.0im
    iz6 = im * z6
    expected6 = -0.5 * (expint(iz6) + expint(-iz6)) + π * im
    @test complex_approx(Ci_complex(z6), expected6)

    # Test case 6: Positive infinity on the real axis
    z7 = Inf + 0.0im
    @test Ci_complex(z7) == 0.0 + 0.0im

    # Test case 7: Negative infinity on the real axis
    z8 = -Inf + 0.0im
    @test Ci_complex(z8) == π * im

    # Test case 8: Complex infinity (negative real and positive imaginary)
    z10 = -Inf + Inf * im
    @test isnan(real(Ci_complex(z10))) && isnan(imag(Ci_complex(z10)))

    # Test case 9: Zero input (real and complex zero)
    z0 = 0.0 + 0.0im
    @test isnan(real(Ci_complex(z0))) && isnan(imag(Ci_complex(z0)))

    z11 = 0.0 + 0.0im
    @test isnan(real(Ci_complex(z11))) && isnan(imag(Ci_complex(z11)))
end

