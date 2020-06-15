using MeshArrays, Statistics, OrdinaryDiffEq

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
   dirIn="run_offflt/"
   prec=Float32
   df=read_flt(dirIn,prec)
end

"""
    example2()

Reproducing `MITgcm/verification/flt_example/` case. This is based on an
extended and modified configuration of the standard MITgcm test case.

```
(df,ref,sol)=example2();

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))
PlotBasic(df,300,100000.0)

using Plots
Plots.plot(sol[1,:],sol[2,:],linewidth=5,lc=:black, title="One Trajectory Example",
xaxis="x",yaxis="y",label="Julia Solution") # legend=false
pl=Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,lc=:red,label="MITgcm Solution")
```
"""
function example2()
   dirIn="flt_example/"
   prec=Float32
   df=read_flt(dirIn,prec)
   uvetc=example2_setup()
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
   ref=ref./uvetc["dx"]
   #
   uInit=[tmp[1,:lon];tmp[1,:lat]]./uvetc["dx"]
   du=fill(0.0,2)
   #
   tspan = (0.0,nSteps*3600.0)
   #prob = ODEProblem(□,uInit,tspan,tmp)
   prob = ODEProblem(⬡,uInit,tspan,uvetc)
   sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
   #
   return df,ref,sol
end

"""
example2_setup()

Read gridded variables from file using MeshArrays and
return result in uvetc Dictionary.
"""
function example2_setup()

   ###### 1) Get gridded variables via MeshArrays.jl


   γ=gcmgrid("flt_example/","PeriodicChannel",1,[(80,42)], [80 42], Float32, read, write)
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

   ## Merge the two dictionaries:

   uvetc=Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1)

   uvetc=merge(uvetc,Γ)

   ## Visualize velocity fields

   mskW=read_mds(γ.path*"hFacW",MeshArray(γ,Float32,nr))
   mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   mskS=read_mds(γ.path*"hFacS",MeshArray(γ,Float32,nr))
   mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)

   msk=Dict("mskW" => mskW, "mskS" => mskS)

   uvetc=merge(uvetc,msk)

end

"""
    example3()

Run simulation over real Ocean domain (-69.5°S to 56.2°N)

```
using MAT
df=example3();

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_pyplot.jl"))
PyPlot.figure(); PlotMapProj(df,3000); gcf()

#include(joinpath(p,"../examples/recipes_Makie.jl"))
#PlotMakie(df,3000,180.)
```
"""
function example3()
   uvetc=example3_setup()

   ii1=0:5:360; ii2=20:2:150;
   n1=length(ii1); n2=length(ii2);
   u0=Array{Float64,2}(undef,(2,n1*n2))
   for i1 in eachindex(ii1); for i2 in eachindex(ii2);
           i=i1+(i2-1)*n1
           u0[1,i]=ii1[i1]
           u0[2,i]=ii2[i2]
   end; end;

   du=fill(0.0,size(u0))
   #⬡(du,u0,uvetc,0.0)
   𝑇 = (0.0,uvetc["t1"]-uvetc["t0"])
   prob = ODEProblem(⬡,u0,𝑇,uvetc)
   sol = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)

   sol[1,:,:]=mod.(sol[1,:,:],360)
   sol[2,:,:]=mod.(sol[2,:,:],180)
   XC=exchange(uvetc["XC"])
   YC=exchange(uvetc["YC"])
   df=postprocess_ODESolution(sol,XC,YC)

   return df
end

"""
example3_setup()

Define gridded variables and return result as Dictionary (`uvetc`).
"""
function example3_setup()

   γ=gcmgrid("llc90_latlon/","PeriodicChannel",1,
                  [(360,178)], [360 178], Float32, read, write)

   Γ=Dict("XC" => read(γ.path*"XC.latlon.data",MeshArray(γ,Float32)),
   "YC" => read(γ.path*"YC.latlon.data",MeshArray(γ,Float32)),
   "XG" => read(γ.path*"XG.data",MeshArray(γ,Float32)),
   "YG" => read(γ.path*"YG.data",MeshArray(γ,Float32)),
   "DXC" => read(γ.path*"DXC.latlon.data",MeshArray(γ,Float32)),
   "DYC" => read(γ.path*"DYC.latlon.data",MeshArray(γ,Float32)) );

   file = matopen(γ.path*"uv_lonlat.mat")
   u=read(file, "u")
   v=read(file, "v")
   close(file)

   u=dropdims(mean(u,dims=3),dims=3)
   v=dropdims(mean(v,dims=3),dims=3)

   u=read(u,MeshArray(γ,Float32))
   v=read(v,MeshArray(γ,Float32));

   u[findall(isnan.(u))]=0.0
   v[findall(isnan.(v))]=0.0

   u0=u; u1=u;
   v0=v; v1=v;

   t0=0.0; t1=86400*366*2.0; dt=3600;

   u0=u0./Γ["DXC"]
   u1=u1./Γ["DXC"]
   v0=v0./Γ["DYC"]
   v1=v1./Γ["DYC"]

   uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;

   uvetc=merge(uvt,Γ);

   nr=50; kk=1;

   mskW=read(γ.path*"hFacW.latlon.data",MeshArray(γ,Float32,nr))
   mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   mskS=read(γ.path*"hFacS.latlon.data",MeshArray(γ,Float32,nr))
   mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)

   msk=Dict("mskW" => mskW, "mskS" => mskS)

   return merge(uvetc,msk)

end
