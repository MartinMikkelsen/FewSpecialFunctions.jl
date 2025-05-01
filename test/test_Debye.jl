@testset "Debye" begin
    @test FewSpecialFunctions.Debye_function(1.0, 0.1) ≈ 0.975278 atol = 1e-3
    @test FewSpecialFunctions.Debye_function(1.0, 0.2) ≈ 0.951111 atol = 1e-4

    @test FewSpecialFunctions.Debye_function(1.0, 1.0) ≈ 0.777505 atol = 1e-4
    @test FewSpecialFunctions.Debye_function(1.0, 1.3) ≈ 0.721173 atol = 1e-5

    @test FewSpecialFunctions.Debye_function(1.0, 2.6) ≈ 0.526375 atol = 1e-5
    @test FewSpecialFunctions.Debye_function(1.0, 2.8) ≈ 0.502682 atol = 1e-5

    @test FewSpecialFunctions.Debye_function(2.0, 1.3) ≈ 0.635800 atol = 1e-5
    @test FewSpecialFunctions.Debye_function(2.0, 2.8) ≈ 0.368324 atol = 1e-5

    @test FewSpecialFunctions.Debye_function(3.0, 0.7) ≈ 0.761859 atol = 1e-5
    @test FewSpecialFunctions.Debye_function(3.0, 2.8) ≈ 0.309995 atol = 1e-5
end