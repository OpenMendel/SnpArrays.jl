using Documenter, SnpArrays

makedocs(
    format = Documenter.HTML(),
    sitename = "SnpArrays.jl",
    authors = "Hua Zhou",
    clean = true,
    debug = true,
    pages = [
        "SnpArrays.jl Tutorial" => "index.md",
        "Linear Algebra Benchmarks" => "linalg.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/SnpArrays.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)

