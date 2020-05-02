module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

using MeshArrays, OrdinaryDiffEq, StatsBase, DataFrames, Random

include("compute.jl")
include("read.jl")
include("examples.jl")
include("update_locations.jl")
include("data_wrangling.jl")

⬡! = IndividualDisplacements.VelComp!
⬡ = IndividualDisplacements.VelComp

export ⬡!, ⬡, ReadDisplacements
export initialize_locations, read_uvetc, postprocess_ODESolution

end # module
