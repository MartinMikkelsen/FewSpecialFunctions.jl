using Documenter, FewSpecialFunctions, Plots

makedocs(sitename="FewSpecialFunctions.jl", pages = ["Home" => "index.md",
"Functions" => "Functions.md"])

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl.git",
)