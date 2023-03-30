using FewSpecialFunctions
using Test

@testset "FewSpecialFunctions.jl" begin

    @test FewSpecialFunctions.regular_Coulomb(0.0, 0.3, 1e-6) ≈ 0.0 atol = 1e-5
    @test FewSpecialFunctions.regular_Coulomb(0.0, -0.3, 1e-6) ≈ 0.0 atol = 1e-5
    @test FewSpecialFunctions.regular_Coulomb(0.0, 2.0, 5.0) ≈ 1.1433 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(1.0, 2.0, 5.0) ≈ 0.99347 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(2.0, 2.0, 5.0) ≈ 0.72125 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(3.0, 2.0, 5.0) ≈ 0.4313 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(1.0, 0.1) ≈ 0.975278 atol = 1e-3
    @test FewSpecialFunctions.Debye_function(1.0, 0.2) ≈ 0.951111 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(1.0, 1.0) ≈ 0.777505 atol = 1e-4
    @test FewSpecialFunctions.Debye_function(1.0, 1.3) ≈ 0.721173 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(1.0, 2.6) ≈ 0.526375 atol = 1e-4
    @test FewSpecialFunctions.Debye_function(1.0, 2.8) ≈ 0.502682 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(2.0, 1.3) ≈ 0.635800 atol = 1e-4
    @test FewSpecialFunctions.Debye_function(2.0, 2.8) ≈ 0.368324 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(3.0, 0.7) ≈ 0.761859 atol = 1e-4
    @test FewSpecialFunctions.Debye_function(3.0, 2.8) ≈ 0.309995 atol = 1e-4

    @test FewSpecialFunctions.Fresnel_S_integral_pi(0.10) ≈ 0.0005236 atol = 1e-4
    @test FewSpecialFunctions.Fresnel_S_integral_pi(0.32) ≈ 0.0171256 atol = 1e-4
    @test FewSpecialFunctions.Fresnel_S_integral_pi(0.70) ≈ 0.1721365 atol = 1e-4

    @test FewSpecialFunctions.Fresnel_C_integral_pi(0.10) ≈ 0.0999975 atol = 1e-4
    @test FewSpecialFunctions.Fresnel_C_integral_pi(0.32) ≈ 0.3191731 atol = 1e-4
    @test FewSpecialFunctions.Fresnel_C_integral_pi(0.70) ≈ 0.6596524 atol = 1e-4

    @test FewSpecialFunctions.Fresnel_S_erf(0.10) ≈ FewSpecialFunctions.Fresnel_S_integral(0.10) atol = 1e-5
    @test FewSpecialFunctions.Fresnel_S_erf(0.32) ≈ FewSpecialFunctions.Fresnel_S_integral(0.32) atol = 1e-5
    @test FewSpecialFunctions.Fresnel_S_erf(0.70) ≈ FewSpecialFunctions.Fresnel_S_integral(0.70) atol = 1e-5

    @test FewSpecialFunctions.Fresnel_C_erf(0.10) ≈ FewSpecialFunctions.Fresnel_C_integral(0.10) atol = 1e-5
    @test FewSpecialFunctions.Fresnel_C_erf(0.32) ≈ FewSpecialFunctions.Fresnel_C_integral(0.32) atol = 1e-5
    @test FewSpecialFunctions.Fresnel_C_erf(0.70) ≈ FewSpecialFunctions.Fresnel_C_integral(0.70) atol = 1e-5


    @test FewSpecialFunctions.Struve(0.0, 4) ≈ 0.1350146 atol = 1e-4
    @test FewSpecialFunctions.Struve(1.0, 4) ≈ 1.0697267 atol = 1e-4
    @test FewSpecialFunctions.Struve(2.0, 4) ≈ 1.2486751 atol = 1e-4

    @test FewSpecialFunctions.Clausen(0) ≈ 0 atol = 1e-4
    @test FewSpecialFunctions.Clausen(π / 3 + 2 * π) ≈ 1.01494160 atol = 1e-4
    @test FewSpecialFunctions.Clausen(-π / 3 + 2 * π) ≈ -1.01494160 atol = 1e-4

    # Test the output for j = 0 and various values of x
    #@test FewSpecialFunctions.FermiDiracIntegral(0, 0) ≈ 0.6931471805599453
    #@test FewSpecialFunctions.FermiDiracIntegral(0, 1) ≈ 1.313261687518223
    #@test FewSpecialFunctions.FermiDiracIntegral(0, 2) ≈ 2.1269280110429727

    # Test the output for j = -1 and various values of x
    #@test FewSpecialFunctions.FermiDiracIntegral(-1, 0) ≈ 0.5
    #@test FewSpecialFunctions.FermiDiracIntegral(-1, 1) ≈ 0.7310585786300049
    #@test FewSpecialFunctions.FermiDiracIntegral(-1, 2) ≈ 0.8807970779778823

    # Test the output for j = 1/2 and x < 1.3
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 0) ≈ 2/(sqrt(π))*0.6780938951530457 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1) ≈ 2/(sqrt(π))*1.3963752806666279 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.2) ≈ 2/(sqrt(π))*1.5863233997463857 atol = 1e-2

    # Test the output for j = 1/2 and x >= 1.3
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 4.1) ≈ 2/(sqrt(π))*5.965800008889902 atol = 1e-1
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 5.2) ≈ 2/(sqrt(π))*8.733390315588004 atol = 1e-1
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 6.6) ≈ 2/(sqrt(π))*11.632113406633252 atol = 1e-1

    # Test the output for j = 3/2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 0.0) ≈ 4/(3*sqrt(π))*1.15280383708879 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1) ≈ 4/(3*sqrt(π))*2.6616826247307124 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.2) ≈ 4/(3*sqrt(π))*3.10869199517456 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.5) ≈ 4/(3*sqrt(π))*0.1075808743944384 atol = 1e-2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.0) ≈ 4/(3*sqrt(π))*0.175800988853926 atol = 1e-2

    #Better tests are needed. Compared to https://npplus.readthedocs.io/en/latest/fermi.html
end
