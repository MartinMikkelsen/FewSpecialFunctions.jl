@testset "Marcum Q-function" begin

    # assume marcumq_test.txt lives in test/data/
    data = open(readdlm, joinpath(@__DIR__, "data", "marcumq_test.txt"))

    @testset "MarcumQ vs MATLAB reference" begin
        @test size(data, 2) == 4

        for r in 1:size(data, 1)
            M, a, b, Q_ref = data[r, :]

            # test the standard MarcumQ(M,a,b)
            @test FewSpecialFunctions.MarcumQ(M, a, b) ≈ Q_ref rtol = 1e-10
        end
    end
    

    # Acta Univ. Sapientiae, Mathematica, 3, 1 (2011) 60–7
    @test FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) ≈ 0.838249985438908 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 0.2, 0.6) ≈ 0.999998670306184 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 0.2, 0.6) ≈  0.999999999927717 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 1.2, 1.6) ≈ 0.501536568390858 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 1.2, 1.6) ≈ 0.994346394491553 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 1.2, 1.6) ≈   0.999944937223540 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9

    @test FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9

    # Compared to marcumq (https://se.mathworks.com/help/signal/ref/marcumq.html)

    @test FewSpecialFunctions.MarcumQ(5, 0.47, 4.85) ≈ 0.0106766402997493 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(5, 1.46, 4.00) ≈ 0.211798804811782  atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 1.27, 4.58) ≈ 0.000931257801666407 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 2.88, 3.28) ≈ 0.773155207859263 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 2.55, 4.67) ≈  0.024112315799424 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 2.27, 3.72) ≈ 0.400088995953665 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(2, 1.97, 0.86) ≈ 0.990345203236692 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 0.10, 1.38) ≈ 0.983869651076909 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(1, 0.29, 4.12) ≈ 0.00028475422733874 atol = 1e-9
    @test FewSpecialFunctions.MarcumQ(4, 0.95, 4.75) ≈ 0.00899422673877906 atol = 1e-9

end
@testset "MarcumQ overloads" begin
    # Test for array inputs
    @testset "Array inputs" begin
        # Test when b is an array
        b_arr = [0.6, 1.6, 2.6]
        q_b_arr = FewSpecialFunctions.MarcumQ(1, 0.2, b_arr)
        @test q_b_arr isa Vector
        @test length(q_b_arr) == length(b_arr)
        @test q_b_arr[1] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test q_b_arr[2] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 1.6) atol = 1e-9
        @test q_b_arr[3] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 2.6) atol = 1e-9
        
        # Test when a is an array
        a_arr = [0.2, 1.2, 2.2]
        q_a_arr = FewSpecialFunctions.MarcumQ(1, a_arr, 0.6)
        @test q_a_arr isa Vector
        @test length(q_a_arr) == length(a_arr)
        @test q_a_arr[1] ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test q_a_arr[2] ≈ FewSpecialFunctions.MarcumQ(1, 1.2, 0.6) atol = 1e-9
        @test q_a_arr[3] ≈ FewSpecialFunctions.MarcumQ(1, 2.2, 0.6) atol = 1e-9
        
        # Test when μ is an array
        μ_arr = [1.0, 5.0, 7.7]
        q_μ_arr = FewSpecialFunctions.MarcumQ(μ_arr, 0.2, 0.6)
        @test q_μ_arr isa Vector
        @test length(q_μ_arr) == length(μ_arr)
        @test q_μ_arr[1] ≈ FewSpecialFunctions.MarcumQ(1.0, 0.2, 0.6) atol = 1e-9
        @test q_μ_arr[2] ≈ FewSpecialFunctions.MarcumQ(5.0, 0.2, 0.6) atol = 1e-9
        @test q_μ_arr[3] ≈ FewSpecialFunctions.MarcumQ(7.7, 0.2, 0.6) atol = 1e-9
    end
    
    # Test for the default μ=1 case
    @testset "Default μ=1" begin
        @test FewSpecialFunctions.MarcumQ(0.2, 0.6) ≈ FewSpecialFunctions.MarcumQ(1, 0.2, 0.6) atol = 1e-9
        @test FewSpecialFunctions.MarcumQ(1.2, 1.6) ≈ FewSpecialFunctions.MarcumQ(1, 1.2, 1.6) atol = 1e-9
        @test FewSpecialFunctions.MarcumQ(2.2, 2.6) ≈ FewSpecialFunctions.MarcumQ(1, 2.2, 2.6) atol = 1e-9
    end
end
