module IndividualDisplacements

greet() = print("Get ready for IndividualDisplacements!")

include("ReadIndDisp.jl")
include("VelocityIndDisp.jl")
include("examples.jl")

export VelComp, VelComp!, VelCopy, ReadGriddedFields, ReadDisplacements

#include("plot_pyplot.jl")
#include("plot_makie.jl")

#export PlotBasic, PlotMapProj, PlotMakie

end # module
