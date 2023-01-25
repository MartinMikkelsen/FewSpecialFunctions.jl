using FewSpecialFunctions
using Test

@testset "FewSpecialFunctions.jl" begin

    @test FewSpecialFunctions.regular_coulomb(0.0,0.3,1e-6) ≈ 0.0 atol=1e-5
    @test FewSpecialFunctions.regular_coulomb(0.0,-0.3,1e-6) ≈ 0.0 atol=1e-5
    @test FewSpecialFunctions.regular_coulomb(0.0,2.0,5.0) ≈ 1.1433 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(1.0,2.0,5.0) ≈ 0.99347 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(2.0,2.0,5.0) ≈ 0.72125 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(3.0,2.0,5.0) ≈ 0.4313 atol=1e-4

    @test FewSpecialFunctions.Debye_function(1.0,0.1) ≈ 0.975278 atol=1e-3
    @test FewSpecialFunctions.Debye_function(1.0,0.2) ≈ 0.951111 atol=1e-4

    @test FewSpecialFunctions.Debye_function(1.0,1.0) ≈ 0.777505 atol=1e-4
    @test FewSpecialFunctions.Debye_function(1.0,1.3) ≈ 0.721173 atol=1e-4

    @test FewSpecialFunctions.Debye_function(1.0,2.6) ≈ 0.526375 atol=1e-4
    @test FewSpecialFunctions.Debye_function(1.0,2.8) ≈ 0.502682 atol=1e-4
 
    @test FewSpecialFunctions.Debye_function(2.0,1.3) ≈ 0.635800 atol=1e-4
    @test FewSpecialFunctions.Debye_function(2.0,2.8) ≈ 0.368324 atol=1e-4

    @test FewSpecialFunctions.Debye_function(3.0,0.7) ≈ 0.761859 atol=1e-4
    @test FewSpecialFunctions.Debye_function(3.0,2.8) ≈ 0.309995 atol=1e-4

    @test FewSpecialFunctions.S_integral_pi(0.10) ≈ 0.0005236 atol=1e-4
    @test FewSpecialFunctions.S_integral_pi(0.32) ≈ 0.0171256 atol=1e-4
    @test FewSpecialFunctions.S_integral_pi(0.70) ≈ 0.1721365 atol=1e-4

    @test FewSpecialFunctions.C_integral_pi(0.10) ≈ 0.0999975 atol=1e-4
    @test FewSpecialFunctions.C_integral_pi(0.32) ≈ 0.3191731 atol=1e-4
    @test FewSpecialFunctions.C_integral_pi(0.70) ≈ 0.6596524 atol=1e-4

    @test FewSpecialFunctions.S_err(0.10) ≈ FewSpecialFunctions.S_integral(0.10) atol=1e-5
    @test FewSpecialFunctions.S_err(0.32) ≈ FewSpecialFunctions.S_integral(0.32) atol=1e-5
    @test FewSpecialFunctions.S_err(0.70) ≈ FewSpecialFunctions.S_integral(0.70) atol=1e-5

    @test FewSpecialFunctions.C_err(0.10) ≈ FewSpecialFunctions.C_integral(0.10) atol=1e-5
    @test FewSpecialFunctions.C_err(0.32) ≈ FewSpecialFunctions.C_integral(0.32) atol=1e-5
    @test FewSpecialFunctions.C_err(0.70) ≈ FewSpecialFunctions.C_integral(0.70) atol=1e-5


    @test FewSpecialFunctions.Struve(0.0,4) ≈ 0.1350146 atol=1e-4
    @test FewSpecialFunctions.Struve(1.0,4) ≈ 1.0697267 atol=1e-4
    @test FewSpecialFunctions.Struve(2.0,4) ≈ 1.2486751 atol=1e-4

    @test FewSpecialFunctions.Struve_alt(0.0,4) ≈ 0.1350146 atol=1e-4
    @test FewSpecialFunctions.Struve_alt(1.0,4) ≈ 1.0697267 atol=1e-4
    @test FewSpecialFunctions.Struve_alt(2.0,4) ≈ 1.2486751 atol=1e-4

    @test FewSpecialFunctions.Clausen(0) ≈ 0 atol=1e-4
    @test FewSpecialFunctions.Clausen(π/3+2*π) ≈ 1.01494160 atol=1e-4
    @test FewSpecialFunctions.Clausen(-π/3+2*π) ≈ -1.01494160 atol=1e-4

end
