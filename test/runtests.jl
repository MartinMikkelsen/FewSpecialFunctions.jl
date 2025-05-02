using FewSpecialFunctions
using Test

@testset "FewSpecialFunctions.jl" begin

    include("test_Clausen.jl")
    include("test_Coulomb.jl")
    include("test_Debye.jl")          
    include("test_FermiDirac.jl")          
    include("test_Fresnel.jl")          
    include("test_hypergeometric.jl")          
    include("test_MarcumQ.jl")          
    include("test_Struve.jl")          
    include("test_parabolic_cylinder.jl")
    
end


