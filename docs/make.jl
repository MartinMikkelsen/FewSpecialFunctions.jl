using Documenter, FewSpecialFunctions

makedocs(
    build = "build",
    sitename = "FewSpecialFunctions.jl",
    pages = ["Home" => "index.md", "Functions" => "Functions.md", "References" => "API.md"],
    format = Documenter.HTML()
)

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "main"
)
