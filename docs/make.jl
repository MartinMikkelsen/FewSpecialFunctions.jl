using Documenter, FewSpecialFunctions

makedocs(source  = "src", build="build", sitename="FewSpecialFunctions.jl", pages = ["Home" => "index.md",
"Functions" => "Functions.md"])

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl.git",
)