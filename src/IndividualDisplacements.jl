module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

using MeshArrays, OrdinaryDiffEq, StatsBase, DataFrames, Random
using NetCDF, Dates, CFTime, CSV

include("compute.jl")
include("read.jl")
include("examples.jl")
include("update_locations.jl")
include("data_wrangling.jl")

⬡! = VelComp!
⬡ = VelComp
□ = VelCopy

export ⬡!, ⬡, □
export initialize_locations, postprocess_ODESolution
export read_flt, read_uvetc, read_drifters

end # module
