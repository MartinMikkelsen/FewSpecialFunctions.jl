using Documenter, FewSpecialFunctions

makedocs(modules = [FewSpecialFunctions], sitename="FewSpecialFunctions.jl", pages = ["Home" => "index.md",
"Functions" => "Functions.md"])

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl.git",
)