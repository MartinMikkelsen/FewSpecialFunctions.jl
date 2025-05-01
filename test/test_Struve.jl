@testset "Struve" begin
    
    @test FewSpecialFunctions.Struve(0.0, 4) ≈ 0.1350146 atol = 1e-4
    @test FewSpecialFunctions.Struve(1.0, 4) ≈ 1.0697267 atol = 1e-4
    @test FewSpecialFunctions.Struve(2.0, 4) ≈ 1.2486751 atol = 1e-4
    
end