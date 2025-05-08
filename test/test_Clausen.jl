using DelimitedFiles

@testset "Clausen" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "Cl1.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test FewSpecialFunctions.Clausen(1,x) ≈ expected rtol=1e-13
        @test FewSpecialFunctions.Clausen(1,-x) ≈ expected rtol=1e-13
        @test FewSpecialFunctions.Clausen(1,x-2pi) ≈ expected rtol=1e-12
        @test FewSpecialFunctions.Clausen(1,x+2pi) ≈ expected rtol=1e-13

    end

    @test FewSpecialFunctions.Clausen(1,0.5) ≈ 0.70358563513784466 rtol=1e-14
    @test FewSpecialFunctions.Clausen(1,1.0) ≈ 0.042019505825368962 rtol=1e-14
    @test FewSpecialFunctions.Clausen(1,2.0) ≈ -0.52054343429085363 rtol=1e-14

end

@testset "f_n function" begin
    # Test with even n (should use sine)
    @test FewSpecialFunctions.f_n(2, 1, π/4) ≈ sin(π/4) rtol=1e-15
    @test FewSpecialFunctions.f_n(2, 2, π/3) ≈ sin(2π/3)/2^2 rtol=1e-15
    @test FewSpecialFunctions.f_n(4, 3, π/6) ≈ sin(3π/6)/3^4 rtol=1e-15
    @test FewSpecialFunctions.f_n(6, 5, 0.5) ≈ sin(5*0.5)/5^6 rtol=1e-15
    
    # Test with odd n (should use cosine)
    @test FewSpecialFunctions.f_n(1, 1, π/4) ≈ cos(π/4) rtol=1e-15
    @test FewSpecialFunctions.f_n(3, 2, π/3) ≈ cos(2π/3)/2^3 rtol=1e-15
    @test FewSpecialFunctions.f_n(5, 3, π/6) ≈ cos(3π/6)/3^5 rtol=1e-15
    @test FewSpecialFunctions.f_n(7, 5, 0.5) ≈ cos(5*0.5)/5^7 rtol=1e-15
    
    # Test special cases
    @test FewSpecialFunctions.f_n(2, 1, 0.0) == 0.0
    @test FewSpecialFunctions.f_n(1, 1, 0.0) == 1.0
    @test FewSpecialFunctions.f_n(2, 1, π) ≈ 0.0 atol=1e-15
    @test FewSpecialFunctions.f_n(1, 1, π) ≈ -1.0 rtol=1e-15
    
    # Test with larger k values
    @test FewSpecialFunctions.f_n(2, 10, 0.1) ≈ sin(10*0.1)/10^2 rtol=1e-14
    @test FewSpecialFunctions.f_n(3, 10, 0.1) ≈ cos(10*0.1)/10^3 rtol=1e-14
end






