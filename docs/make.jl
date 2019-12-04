using Documenter
using IndividualDisplacements

makedocs(
    sitename = "IndividualDisplacements",
    format = Documenter.HTML(),
    modules = [IndividualDisplacements]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gaelforget/IndividualDisplacements.jl.git",
)
