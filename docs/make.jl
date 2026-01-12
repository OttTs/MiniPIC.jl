using Documenter, MiniPIC

makedocs(
    modules = [MiniPIC],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "MiniPIC.jl",
    authors  = "Tobias Ott",
    pages = [
        "Home" => "index.md",
        "Physical Model" => "physical-model.md",
        "Numerics" => "numerics.md",
        "User guide" => "userguide.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo   = "github.com/OttTs/MiniPIC.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true
)