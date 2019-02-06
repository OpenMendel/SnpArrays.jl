using Documenter, SnpArrays

makedocs(
    format = :html,
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
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)

