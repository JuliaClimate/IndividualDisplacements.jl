
# Download grid, velocity, etc files into `examples/`

using OceanStateEstimation, IndividualDisplacements

p=dirname(pathof(IndividualDisplacements))
cd(joinpath(p,"../examples/"))
include(joinpath(p,"../examples/helper_functions.jl"))
get_grid_if_needed()

q=dirname(pathof(OceanStateEstimation))
run(`cp $q/../examples/nctiles_climatology.csv $p/../examples/`)
run(`mkdir nctiles_climatology`)
get_from_dataverse("UVELMASS","nctiles_climatology/")
get_from_dataverse("VVELMASS","nctiles_climatology/")
