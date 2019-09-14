using Documenter
using IndividualDisplacement

makedocs(
    sitename = "IndividualDisplacement",
    format = :html,
    modules = [IndividualDisplacement]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
