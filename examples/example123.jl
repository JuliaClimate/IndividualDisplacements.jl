using IndividualDisplacements, MeshArrays, MITgcmTools, OrdinaryDiffEq, Statistics

"""
    example1()

Pre-computed global ocean case. Here we just re-read data from a file produced
earlier, rather than computing trajectories as in the other examples.

```
df=example1()

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_pyplot.jl"))
PyPlot.figure(); PlotMapProj(df,300); gcf()
```
"""
function example1()
   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/run_offflt/")
   prec=Float32
   df=read_flt(dirIn,prec)
end

"""
    example2()

Reproducing `MITgcm/verification/flt_example/` case. This is based on an
extended and modified configuration of the standard MITgcm test case.

```
(𝐼,df,ref)=example2();

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))
PlotBasic(df,300,100000.0)

using Plots
Plots.plot(𝐼.🔴.x,𝐼.🔴.y,linewidth=5,lc=:black, title="One Trajectory Example",
xaxis="x",yaxis="y",label="Julia Solution") # legend=false
pl=Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,lc=:red,label="MITgcm Solution")
```
"""
function example2()
   𝑃,Γ=example2_setup()
   (xy,df,ref,nSteps)=example2_xy(𝑃)

   𝑃.𝑇[:] = [0.0,nSteps*3600.0]
   solv(prob) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

   tr = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
   
   𝐼 = Individuals{Float64}(📌=xy[:,:], 🔴=tr, 🚄 = dxy_dt, ∫ = solv, 🔧 = postprocess_xy, 𝑃=𝑃)
   𝑇=(0.0,𝐼.𝑃.𝑇[2])
   ∫!(𝐼,𝑇)

   return 𝐼, df,ref
end

"""
example2_xy()

Read MITgcm/pkg/flt output
"""
function example2_xy(𝑃)
dirIn=IndividualDisplacements.flt_example
prec=Float32
df=read_flt(dirIn,prec)
#
tmp=df[df.ID .== 200, :]
nSteps=Int32(tmp[end,:time]/3600)-2
ref=transpose([tmp[1:nSteps,:lon] tmp[1:nSteps,:lat]])
dx=5000.0
maxLon=80*dx
maxLat=42*dx
for i=1:nSteps-1
    ref[1,i+1]-ref[1,i]>maxLon/2 ? ref[1,i+1:end]-=fill(maxLon,(nSteps-i)) : nothing
    ref[1,i+1]-ref[1,i]<-maxLon/2 ? ref[1,i+1:end]+=fill(maxLon,(nSteps-i)) : nothing
    ref[2,i+1]-ref[2,i]>maxLat/2 ? ref[2,i+1:end]-=fill(maxLat,(nSteps-i)) : nothing
    ref[2,i+1]-ref[2,i]<-maxLat/2 ? ref[2,i+1:end]+=fill(maxLat,(nSteps-i)) : nothing
end
ref=ref./dx
xy=[tmp[1,:lon];tmp[1,:lat]]./dx
return xy,df,ref,nSteps
end

"""
example2_setup()

Read gridded variables from file using MeshArrays and
return result in uvetc Dictionary.
"""
function example2_setup()

   ###### 1) Get gridded variables via MeshArrays.jl

   dirIn=IndividualDisplacements.flt_example
   γ=gcmgrid(dirIn,"PeriodicChannel",1,[(80,42)], [80 42], Float32, read, write)
   nr=8

   ## Put grid variables in a dictionary:

   Γ=Dict("XC" => read_mds(γ.path*"XC",MeshArray(γ,Float32)),
   "YC" => read_mds(γ.path*"YC",MeshArray(γ,Float32)),
   "XG" => read_mds(γ.path*"XG",MeshArray(γ,Float32)),
   "YG" => read_mds(γ.path*"YG",MeshArray(γ,Float32)),
   "dx" => 5000.0)

   ## Put velocity fields in a dictionary:

   t0=0.0 #approximation / simplification
   t1=18001.0*3600.0

   u0=read_mds(γ.path*"U.0000000001",MeshArray(γ,Float32,nr))
   u1=read_mds(γ.path*"U.0000018001",MeshArray(γ,Float32,nr))
   v0=read_mds(γ.path*"V.0000000001",MeshArray(γ,Float32,nr))
   v1=read_mds(γ.path*"V.0000018001",MeshArray(γ,Float32,nr))

   kk=3 #3 to match -1406.25 in pkg/flt output
   u0=u0[:,kk]; u1=u1[:,kk];
   v0=v0[:,kk]; v1=v1[:,kk];

   u0=u0./Γ["dx"]
   u1=u1./Γ["dx"]
   v0=v0./Γ["dx"]
   v1=v1./Γ["dx"]

   ## Visualize velocity fields

   mskW=read_mds(γ.path*"hFacW",MeshArray(γ,Float32,nr))
   mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   mskS=read_mds(γ.path*"hFacS",MeshArray(γ,Float32,nr))
   mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)
   Γ=merge(Γ,Dict("mskW" => mskW, "mskS" => mskS))

   𝑃=𝑃_Array2D{eltype(u0)}(u0[1], u1[1], v0[1], v1[1], [t0,t1])
   return 𝑃,Γ
end
