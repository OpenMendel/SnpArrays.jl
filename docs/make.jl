using Documenter, SnpArrays

run(`jupyter nbconvert --to markdown ./docs/snparray.ipynb --output ./docs/src/man/snparray.md`)
makedocs()
deploydocs(
  deps   = Deps.pip("mkdocs", "python-markdown-math", "jupyter"),
  repo   = "github.com:OpenMendel/SnpArrays.jl.git",
  julia  = "0.4",
  osname = "osx"
  )
