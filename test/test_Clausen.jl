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