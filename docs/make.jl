using Documenter, FewSpecialFunctions

makedocs(
    build = "build",
    sitename = "FewSpecialFunctions.jl", 
    pages = ["Home" => "index.md", "Functions" => "Functions.md", "API" => "API.md"],
    format = Documenter.HTML()
)


deploydocs(
    repo = "https://github.com/MartinMikkelsen/FewSpecialFunctions.jl",
    target = "build",       
    branch = "gh-pages",
    devbranch = "main"
)