using Documenter, FewSpecialFunctions

makedocs(sitename="FewSpecialFunctions.jl", doctest = true, pages=["Home" => "index.md",
"Functions" => "Functions.md"])

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl.git",
)