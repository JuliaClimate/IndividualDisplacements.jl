module IndividualDisplacement

greet() = print("Get ready for IndividualDisplacement!")

include("ReadIndDisp.jl")
include("PlotIndDisp.jl")
include("VelocityIndDisp.jl")
#export ReadGriddedFields, ReadDisplacements
#export PlotBasic, PlotMapProj

# calling sequences

if false
   #for global ocean case, e.g.:
   dirIn="run_offflt/"
   prec=Float32
   df=IndividualDisplacement.ReadDisplacements(dirIn,prec)
   PyPlot.figure()
   IndividualDisplacement.PlotMapProj(df,300)
end

if false
   #for flt_example case, e.g.:
   dirIn="run.long2/"
   prec=Float32
   df=IndividualDisplacement.ReadDisplacements(dirIn,prec)
   PyPlot.figure()
   IndividualDisplacement.PlotBasic(df,300)
   gcf()
   #
   tmp=df[df.ID .== 200, :]
   nSteps=Int32(tmp[end,:time]/3600)-2
   ref=transpose([tmp[1:nSteps,:lon] tmp[1:nSteps,:lat]])
   maxLon=80*5.e3
   maxLat=42*5.e3
   for i=1:nSteps-1
       ref[1,i+1]-ref[1,i]>maxLon/2 ? ref[1,i+1:end]-=fill(maxLon,(nSteps-i)) : nothing
       ref[1,i+1]-ref[1,i]<-maxLon/2 ? ref[1,i+1:end]+=fill(maxLon,(nSteps-i)) : nothing
       ref[2,i+1]-ref[2,i]>maxLat/2 ? ref[2,i+1:end]-=fill(maxLat,(nSteps-i)) : nothing
       ref[2,i+1]-ref[2,i]<-maxLat/2 ? ref[2,i+1:end]+=fill(maxLat,(nSteps-i)) : nothing
   end
   #
   comp_vel=IndividualDisplacement.VelComp
   get_vel=IndividualDisplacement.VelCopy
   uInit=[tmp[1,:lon];tmp[1,:lat]]
   du=fill(0.0,2)
   #
   uvetc=IndividualDisplacement.ReadGriddedFields()
   tspan = (0.0,nSteps*3600.0)
   #prob = ODEProblem(get_vel,uInit,tspan,tmp)
   prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
   sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
   #
   Plots.plot(sol[1,:],sol[2,:],linewidth=5,title="Using Recomputed Velocities",
   xaxis="lon",yaxis="lat",label="Julia Solution") # legend=false
   Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,label="MITgcm Solution")
end

end # module
