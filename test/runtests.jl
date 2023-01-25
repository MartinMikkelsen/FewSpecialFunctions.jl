using FewSpecialFunctions
using Test, Base

@testset "FewSpecialFunctions.jl" begin
    @test FewSpecialFunctions.regular_coulomb(0.0,0.3,1e-6) ≈ 0.0 atol=1e-5
    @test FewSpecialFunctions.regular_coulomb(0.0,-0.3,1e-6) ≈ 0.0 atol=1e-5
    @test FewSpecialFunctions.regular_coulomb(0.0,2.0,5.0) ≈ 1.1433 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(1.0,2.0,5.0) ≈ 0.99347 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(2.0,2.0,5.0) ≈ 0.72125 atol=1e-4
    @test FewSpecialFunctions.regular_coulomb(3.0,2.0,5.0) ≈ 0.4313 atol=1e-4
    #@test FewSpecialFunctions.Debye_function(1.0,0) ≈ 1.0 atol=1e-5
end
