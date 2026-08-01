using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        FewSpecialFunctions; persistent_tasks = false
    )
end
