using Documenter, BEDFiles

makedocs(
    format = :html,
    sitename = "BEDFiles"
)

deploydocs(
    repo   = "github.com/dmbates/BEDFiles.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)

