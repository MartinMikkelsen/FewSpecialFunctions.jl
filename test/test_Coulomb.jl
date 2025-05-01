@testset "Coulomb" begin
    
    @test FewSpecialFunctions.regular_Coulomb(0.0, 0.3, 1e-6) ≈ 0.0 atol = 1e-5
    @test FewSpecialFunctions.regular_Coulomb(0.0, -0.3, 1e-6) ≈ 0.0 atol = 1e-5
    @test FewSpecialFunctions.regular_Coulomb(0.0, 2.0, 5.0) ≈ 1.1433 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(1.0, 2.0, 5.0) ≈ 0.99347 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(2.0, 2.0, 5.0) ≈ 0.72125 atol = 1e-4
    @test FewSpecialFunctions.regular_Coulomb(3.0, 2.0, 5.0) ≈ 0.4313 atol = 1e-4
    
end