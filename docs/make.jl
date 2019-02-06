using Documenter, SnpArrays

makedocs(
    format = Documenter.HTML(),
    sitename = "SnpArrays.jl",
    authors = "Hua Zhou",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/SnpArrays.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)

