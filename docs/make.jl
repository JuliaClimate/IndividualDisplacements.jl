using Documenter, Literate
using IndividualDisplacements

#download data dependencies if needed
IndividualDisplacements.get_ecco_velocity_if_needed();
IndividualDisplacements.get_occa_velocity_if_needed();

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/")
notebooks = joinpath(src, "notebooks")

execute = true # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks

lst1 = ["solid_body_rotation","random_flow_field","global_ocean_circulation","three_dimensional_ocean","detailed_look","particle_cloud"]
#lst2 = ["solid_body_rotation","random_flow_field"]
lst2 = ["none"]
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
i=findall((occursin).("solid_body_rotation",p)); 
j=findall((occursin).("random_flow_field",p)); 
p_tu=[p[j];p[i]]
i=findall((occursin).("detailed_look",p)); 
j=findall((occursin).("particle_cloud",p)); 
p_mi=[p[i];p[j]]

makedocs(
    sitename = "IndividualDisplacements",
    format = Documenter.HTML(),
    pages = [
		"Introduction" => "index.md",
        "User Guide" => "workflow.md",
		"Tool Box" => "API.md",
        "Example Guide" => "examples.md",
        "Tutorial Examples" => p_tu, 
		"Real Ocean Cases" => pages("worldwide"),
        "MITgcm Examples" => p_mi], 
    modules = [IndividualDisplacements]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaClimate/IndividualDisplacements.jl.git",
)
