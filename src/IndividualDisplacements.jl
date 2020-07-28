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

export ⬡!, ⬡, □, duvw
export initialize_gridded, initialize_lonlat, randn_lonlat
export postprocess_lonlat, postprocess_xy
export read_flt, read_mds, read_uvetc, read_drifters

end # module
