using DelimitedFiles
using Test

@testset "Debye" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "debye_test.txt"))

    @testset "Debye_function vs DebyeFunctions.jl" begin
        for r in 1:size(data, 1)
            x, n, Q_ref = data[r, :]
            @test FewSpecialFunctions.debye_function(n, 1.0, x) ≈ Q_ref atol= 1e-4
        end
    end

    @test FewSpecialFunctions.debye_function(2.0, 1.0, 5.0) ≈ 0.172329034857624782145 atol=1e-6
    @test FewSpecialFunctions.debye_function(5.0, 1.0, 4.5) ≈ 0.10164118339698890968 atol=1e-15
    @test FewSpecialFunctions.debye_function(7.0, 1.0, 0.8) ≈ 0.69112406526865230673 atol=1e-14
    @test FewSpecialFunctions.debye_function(9.0, 1.0, 3.4) ≈ 0.15413773867789254146 atol=1e-15
    @test FewSpecialFunctions.debye_function(12.0, 1.0, 5.4) ≈ 0.03618849233828133 atol=1e-15
    @test FewSpecialFunctions.debye_function(15.0, 1.0, 2.4) ≈ 0.2661409156647294951955 atol=1e-15
    @test FewSpecialFunctions.debye_function(20.0, 1.0, 1.24) ≈ 0.523361585088859680745 atol=1e-11
    @test FewSpecialFunctions.debye_function(25.0, 1.0, 4.2) ≈ 0.07296496706218587 atol=1e-15
    @test FewSpecialFunctions.debye_function(30.0, 1.0, 3.42) ≈ 0.1258426106590655660781 atol=1e-13

end
