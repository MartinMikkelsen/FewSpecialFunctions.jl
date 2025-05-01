@testset "Clausen" begin
    
    @test FewSpecialFunctions.Clausen(0) ≈ 0 atol = 1e-4
    @test FewSpecialFunctions.Clausen(π / 3 + 2 * π) ≈ 1.01494160 atol = 1e-4
    @test FewSpecialFunctions.Clausen(-π / 3 + 2 * π) ≈ -1.01494160 atol = 1e-4

end