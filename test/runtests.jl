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

    # Comparing to https://npplus.readthedocs.io/en/latest/fermi.html
    # Test the output for j = -1/2
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 0.0) ≈ 1.0721549299400754 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.0) ≈ 1.8204113571471041 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.2) ≈ 1.9785617633438695 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 4.1) ≈ 3.928454737099184 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 5.2) ≈ 4.477432715418454 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 6.6) ≈ 5.082787981164429 atol = 1e-12

    # Test the output for j = 1/2 
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 0.0) ≈ 0.6780938951530457 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.0) ≈ 1.3963752806666279 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.2) ≈ 1.5863233997463857 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 4.1) ≈ 5.965800008889902 atol = 1e-5
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 5.2) ≈ 8.28102917922544 atol = 1e-6
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 6.6) ≈ 11.632113406633252 atol = 1e-7

    # Test the output for j = 3/2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 0.0) ≈ 1.15280383708879 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247307124 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.2) ≈ 3.10869199517456 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.5) ≈ 0.1075808743944384 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.0) ≈ 0.175800988853926 atol = 1e-12

    # Test the output for j = 5/2
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 0.0) ≈ 3.0825860828379246 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.0) ≈ 7.626535355004442 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.2) ≈ 9.066754659807005 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.5) ≈ 0.27085618652992494 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.0) ≈ 0.4445544534586879 atol = 1e-13

    # Comparing to https://github.com/scott-maddox/fdint

    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 0.0) ≈ 1.0721549299401913 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.0) ≈ 1.8204113571469622 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.2) ≈ 1.9785617633438086 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 4.1) ≈ 3.928454737096438 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 5.2) ≈ 4.477432715417047 atol = 1e-12
    # @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 6.6) ≈ 5.082787981168395 atol = 1e-12

    # Test the output for j = 1/2 
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 0.0) ≈ 0.6780938951531009 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.0) ≈ 1.396375280666564 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.2) ≈ 1.586323399746276 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 4.1) ≈ 5.965800008886879 atol = 1e-5
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 5.2) ≈ 8.281029179226932 atol = 1e-6
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 6.6) ≈ 11.632113406638712 atol = 1e-5

    # Test the output for j = 3/2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 0.0) ≈ 1.1528038370883613 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247320042 atol = 1e-10
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.2) ≈ 3.108691995173995 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.5) ≈ 0.10758087439445904 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.0) ≈ 0.1758009888540129 atol = 1e-12

    # Test the output for j = 5/2
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 0.0) ≈ 3.0825860828374183 atol = 1e-12
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.0) ≈ 7.626535355005596 atol = 1e-10
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.2) ≈ 9.066754659806907 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.5) ≈ 0.2708561865299367 atol = 1e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.0) ≈ 0.44455445345876304 atol = 1e-13

    end


