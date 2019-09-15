module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

include("ReadIndDisp.jl")
include("PlotIndDisp.jl")
include("VelocityIndDisp.jl")

export ReadGriddedFields, ReadDisplacements
export PlotBasic, PlotMapProj

end # module
