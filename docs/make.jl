using Documenter, SnpArrays

run(`jupyter nbconvert --to markdown snparray.ipynb --output ./src/man/snparray.md`)
makedocs()
deploydocs(
  deps   = Deps.pip("mkdocs", "python-markdown-math", "jupyter"),
  repo   = "github.com:OpenMendel/SnpArrays.jl.git",
  julia  = "0.4",
  osname = "osx"
  )
