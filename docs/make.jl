using Documenter, SnpArrays

makedocs()
deploydocs(
  deps   = Deps.pip("mkdocs", "python-markdown-math"),
  repo   = "github.com:OpenMendel/SnpArrays.jl.git",
  julia  = "0.4",
  osname = "osx"
  )
