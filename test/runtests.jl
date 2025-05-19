using FewSpecialFunctions
using Test
using DelimitedFiles  

@testset "FewSpecialFunctions.jl" begin

    include("Aqua.jl")
    include("test_Coulomb.jl")
    include("test_Clausen.jl")
    include("test_Debye.jl")          
    include("test_FermiDirac.jl")          
    include("test_Fresnel.jl")          
    include("test_MarcumQ.jl")          
    include("test_Struve.jl")          
    include("test_parabolic_cylinder.jl")

end


