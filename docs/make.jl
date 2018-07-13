using Documenter, SnpArrays

ENV["DOCUMENTER_DEBUG"] = "true"
makedocs()
deploydocs(
  deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
  repo   = "github.com:OpenMendel/SnpArrays.jl.git",
  julia  = "0.6",
  osname = "linux"
  )
