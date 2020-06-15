module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

using MeshArrays, OrdinaryDiffEq, StatsBase, DataFrames, Random
using NetCDF, Dates, CFTime, CSV

include("compute.jl")
include("read.jl")
include("update_locations.jl")
include("data_wrangling.jl")

⬡! = VelComp!
⬡ = VelComp
□ = VelCopy

export ⬡!, ⬡, □
export initialize_grid_locations, initialize_random_locations
export setup_periodic_domain, randn_lonlat
export postprocess_ODESolution, postprocess_ODESolution_simple
export read_flt, read_mds, read_uvetc, read_drifters

end # module
