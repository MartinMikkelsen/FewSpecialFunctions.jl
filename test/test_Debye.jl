using DelimitedFiles
using Test

@testset "Debye" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "debye_test.txt"))

    @testset "Debye_function vs DebyeFunctions.jl" begin
        for r in 1:size(data, 1)
            x, n, Q_ref = data[r, :]
            @test FewSpecialFunctions.debye_function(n, 1.0, x) ≈ Q_ref atol = 1.0e-4
        end
    end

    @test FewSpecialFunctions.debye_function(2.0, 1.0, 5.0) ≈ 0.172329034857624782145 atol = 1.0e-6
    @test FewSpecialFunctions.debye_function(5.0, 1.0, 4.5) ≈ 0.10164118339698890968 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(7.0, 1.0, 0.8) ≈ 0.69112406526865230673 atol = 1.0e-14
    @test FewSpecialFunctions.debye_function(9.0, 1.0, 3.4) ≈ 0.15413773867789254146 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(12.0, 1.0, 5.4) ≈ 0.03618849233828133 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(15.0, 1.0, 2.4) ≈ 0.2661409156647294951955 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(20.0, 1.0, 1.24) ≈ 0.523361585088859680745 atol = 1.0e-11
    @test FewSpecialFunctions.debye_function(25.0, 1.0, 4.2) ≈ 0.07296496706218587 atol = 1.0e-15
    @test FewSpecialFunctions.debye_function(30.0, 1.0, 3.42) ≈ 0.1258426106590655660781 atol = 1.0e-13

end

@testset "Debye function method variants" begin
    # Test real number conversion
    @test FewSpecialFunctions.debye_function(2, 1, 5) ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 5.0)
    @test FewSpecialFunctions.debye_function(Int8(3), Float32(1.0), 2) ≈ FewSpecialFunctions.debye_function(3.0, 1.0, 2.0)

    # Test array input for x
    x_array = [1.0, 2.0, 3.0]
    result = FewSpecialFunctions.debye_function(2.0, 1.0, x_array)
    @test length(result) == length(x_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 1.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 2.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)

    # Test array input for β
    β_array = [1.0, 1.5, 2.0]
    result = FewSpecialFunctions.debye_function(2.0, β_array, 3.0)
    @test length(result) == length(β_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.5, 3.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(2.0, 2.0, 3.0)

    # Test array input for n
    n_array = [1.0, 2.0, 3.0]
    result = FewSpecialFunctions.debye_function(n_array, 1.0, 3.0)
    @test length(result) == length(n_array)
    @test result[1] ≈ FewSpecialFunctions.debye_function(1.0, 1.0, 3.0)
    @test result[2] ≈ FewSpecialFunctions.debye_function(2.0, 1.0, 3.0)
    @test result[3] ≈ FewSpecialFunctions.debye_function(3.0, 1.0, 3.0)

    # Test default n=1.0
    @test FewSpecialFunctions.debye_function(1.5, 2.0) ≈ FewSpecialFunctions.debye_function(1.0, 1.5, 2.0)
end
