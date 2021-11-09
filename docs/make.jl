using Documenter, Literate, PlutoSliderServer
using IndividualDisplacements, OceanStateEstimation
import CairoMakie as Mkie

#download data dependencies if needed
IndividualDisplacements.flt_example_download()
OceanStateEstimation.get_ecco_velocity_if_needed();
OceanStateEstimation.get_occa_velocity_if_needed();

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/jupyter/")
notebooks = joinpath(src, "notebooks")

execute = true # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks

lst1 = ["detailed_look","particle_cloud","global_ocean_circulation","three_dimensional_ocean"]
lst2 = ["detailed_look","particle_cloud","global_ocean_circulation","three_dimensional_ocean"]
tst1(x) = !isempty(lst1) && Bool(sum(isequal.(x, lst1)))
tst2(x) = !isempty(lst2) && Bool(sum(isequal.(x, lst2)))

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
		"Introduction" => "index.md",
        "User Guide" => "workflow.md",
        "Examples" => "examples.md", 
		"Tool Box" => "API.md"],
        doctest = false,
    modules = [IndividualDisplacements]
)

pth_in = joinpath(@__DIR__, "..","examples")
pth_out = joinpath(@__DIR__, "build","examples")
lst=("solid_body_rotation.jl","random_flow_field.jl","global_ocean_circulation.jl","three_dimensional_ocean.jl")
subpth=("basics","basics","worldwide","worldwide")
for ii in 1:length(lst)
    fil_in=joinpath(pth_in,subpth[ii],lst[ii])
    fil_out=joinpath(pth_out,lst[ii][1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    #cp(fil_in[1:end-2]*"html",fil_out)
end

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaClimate/IndividualDisplacements.jl.git",
)
