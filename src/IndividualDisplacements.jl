module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

using MeshArrays, OrdinaryDiffEq, StatsBase, DataFrames, Random

include("compute.jl")
include("read.jl")
include("examples.jl")
include("update_locations.jl")
include("data_wrangling.jl")

⬡! = VelComp!
⬡ = VelComp
□ = VelCopy

export ⬡!, ⬡, □, ReadDisplacements, read_uvetc
export initialize_locations, postprocess_ODESolution

end # module
