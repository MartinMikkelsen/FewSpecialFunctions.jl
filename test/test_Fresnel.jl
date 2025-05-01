@testset "Fresnel" begin
    
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

end