using Documenter, Literate
using IndividualDisplacements, MITgcmTools, OceanStateEstimation

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/")
notebooks = joinpath(src, "notebooks")

execute = true # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks

lst1 = ["solid_body_rotation","random_flow_field","global_ocean_circulation","detailed_look","particle_cloud"]
lst2 = ["solid_body_rotation","random_flow_field","global_ocean_circulation","detailed_look","particle_cloud"]
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

p=pages("basics"); np=length(p)
i=findall((!occursin).("detailed_look",p)); i=[i;setdiff(1:np,i)]; p=p[i]
i=findall((occursin).("solid_body_rotation",p)); i=[i;setdiff(1:np,i)]; p=p[i]
p_BE=p


makedocs(
    sitename = "IndividualDisplacements",
    format = Documenter.HTML(),
    pages = [
		"Introduction" => "index.md",
                "API" => "workflow.md",
		"Examples" => "examples.md",
		"Real Ocean" => pages("worldwide"),
                "Other Examples" => p_BE,
		"Tool Boxes" => "API.md"],
    modules = [IndividualDisplacements]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaClimate/IndividualDisplacements.jl.git",
)
