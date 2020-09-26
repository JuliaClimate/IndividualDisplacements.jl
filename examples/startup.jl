# Download grid, velocity, etc files into `examples/`

using OceanStateEstimation, IndividualDisplacements

p=dirname(pathof(IndividualDisplacements))
cd(joinpath(p,"../examples/"))
include(joinpath(p,"../examples/helper_functions.jl"))
get_ll360_grid_if_needed(); get_occa_velocity_if_needed();
get_llc90_grid_if_needed(); get_ecco_velocity_if_needed();

