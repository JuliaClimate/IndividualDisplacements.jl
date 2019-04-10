module IndividualDisplacement

greet() = print("Get ready for IndividualDisplacement!")

include("ReadIndDisp.jl")
include("PlotIndDisp.jl")
include("VelocityIndDisp.jl")
#export ReadGriddedFields, ReadDisplacements
#export PlotBasic, PlotMapProj

# example of calling sequence

if false
   #for flt_example case, e.g.:
   dirIn="run.long2/"
   prec=Float32
   df=IndividualDisplacement.ReadDisplacements(dirIn,prec);
   PyPlot.figure()
   IndividualDisplacement.PlotBasic(df,300);
   #
   tmp=df[df.ID .== 200, :]
   tmp[1:4,:]

end

if false
   #for global ocean case, e.g.:
   dirIn="run_offflt/"
   prec=Float32
   df=IndividualDisplacement.ReadDisplacements(dirIn,prec)
   PyPlot.figure()
   IndividualDisplacement.PlotMapProj(df,300)
end

end # module
