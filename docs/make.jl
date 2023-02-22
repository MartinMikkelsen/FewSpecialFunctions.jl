using Documenter, FewSpecialFunctions

makedocs(
    source  = "src", 
    workdir = "build", 
    sitename = "FewSpecialFunctions.jl", 
    pages = ["Home" => "index.md", "Functions" => "Functions.md", "API" => "API.md"],
    format = Documenter.HTML()
)

deploydocs(
    repo = "https://github.com/MartinMikkelsen/FewSpecialFunctions.jl",
)
