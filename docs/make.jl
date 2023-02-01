using Documenter, FewSpecialFunctions

makedocs(sitename="FewSpecialFunctions.jl", pages = ["Home" => "index.md",modules = [FewSpecialFunctions],
"Functions" => "Functions.md"])

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl.git",
)