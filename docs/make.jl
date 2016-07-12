using Documenter, SnpArrays, IJulia
import Conda

if success(`which jupyter`)
  jupyter_program = "jupyter"
else
  jupyter_program = joinpath(Conda.SCRIPTDIR, "jupyter")
end
run(`$jupyter_program nbconvert --to markdown ./docs/snparray.ipynb --output ./docs/src/man/snparray.md`)
makedocs()
deploydocs(
  deps   = Deps.pip("mkdocs", "python-markdown-math"),
  repo   = "github.com:OpenMendel/SnpArrays.jl.git",
  julia  = "0.4",
  osname = "osx"
  )
