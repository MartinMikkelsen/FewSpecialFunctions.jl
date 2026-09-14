using Documenter, FewSpecialFunctions

DocMeta.setdocmeta!(FewSpecialFunctions, :DocTestSetup, :(using FewSpecialFunctions); recursive = true)

makedocs(
    build = "build",
    sitename = "FewSpecialFunctions.jl",
    modules = [FewSpecialFunctions],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Functions" => "Functions.md",
        "References" => "API.md",
    ],
    format = Documenter.HTML()
)

deploydocs(
    repo = "github.com/MartinMikkelsen/FewSpecialFunctions.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "main"
)
