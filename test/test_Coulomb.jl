@testset "Coulomb" begin
    
@testset "CoulombF" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "CoulombF.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        ℓ        = row[1]
        η        = row[2]
        ρ        = row[3]
        expected = row[4]

        @test real.(FewSpecialFunctions.F(ℓ,η,ρ)) ≈ expected atol=1e-8

    end

end

@testset "CoulombG" begin

    data = open(readdlm, joinpath(@__DIR__, "data", "CoulombG.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        ℓ        = row[1]
        η        = row[2]
        ρ        = row[3]
        expected = row[4]

        @test real.(FewSpecialFunctions.G(ℓ,η,ρ)) ≈ expected atol=1e-5

    end

end
    
end