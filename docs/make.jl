using Documenter, Literate
using IndividualDisplacements

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/")
notebooks = joinpath(src, "notebooks")

execute = true # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks

lst1 = ["detailed_look","particle_cloud","solid_body_rotation","random_flow_field","global_ocean_circulation"]
lst2 = ["detailed_look","particle_cloud","solid_body_rotation","random_flow_field"]
tst1(x) = Bool(sum(isequal.(x, lst1)))
tst2(x) = Bool(sum(isequal.(x, lst2)))

for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
	tst1(splitext(file)[1]) || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    nb && Literate.notebook(ipath, notebooks, execute = execute*tst2(splitext(file)[1]))
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) = [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

makedocs(
    sitename = "IndividualDisplacements",
    format = Documenter.HTML(),
    pages = [
		"Home" => "index.md",
		"List Of Examples" => "examples.md",
		"Basic Examples" => pages("basics"),
		"Global Examples" => pages("worldwide"),
		"API Guide" => "API.md"],
    modules = [IndividualDisplacements]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaClimate/IndividualDisplacements.jl.git",
)
