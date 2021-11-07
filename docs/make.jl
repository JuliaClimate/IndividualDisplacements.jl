using Documenter, Literate, PlutoSliderServer
using IndividualDisplacements, OceanStateEstimation

#download data dependencies if needed
IndividualDisplacements.flt_example_download()
OceanStateEstimation.get_ecco_velocity_if_needed();
OceanStateEstimation.get_occa_velocity_if_needed();

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/")
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

p=pages("basics"); np=length(p)
i=findall((occursin).("detailed_look",p)); 
j=findall((occursin).("particle_cloud",p)); 
p_mi=[p[i];p[j]]

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

pth = joinpath(@__DIR__, "build","examples")
lst=("solid_body_rotation.jl","random_flow_field.jl")
for i in lst
    fil_in=joinpath(@__DIR__,"..","examples","basics",i)
    fil_out=joinpath(pth,i[1:end-2]*"html")
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
