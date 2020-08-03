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
   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/flt_example/")
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

   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/flt_example/")
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
using MAT, NetCDF
df=example3("OCCA")

p=dirname(pathof(IndividualDisplacements))

include(joinpath(p,"../examples/recipes_plots.jl"))
PlotBasic(df,1000,90.0)

#include(joinpath(p,"../examples/recipes_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,3000); gcf()

#include(joinpath(p,"../examples/recipes_Makie.jl"))
#PlotMakie(df,3000,180.)

n=maximum(df.ID)
nt=size(df,1)/n
dt=maximum(df.t)/(nt-1)
@gif for t in 0:nt-1
   scatter_zcolor(df,t*dt,df.z,(0,10))
end
```
"""
function example3(nam::String="OCCA";backward_in_time::Bool=false,
   lon_rng=(-165.0,-145.0), lat_rng=(25.0,35.0), z_init=1.0)
   if nam=="OCCA"
      uvetc=OCCA_setup(backward_in_time)
   elseif nam=="LLC90"
      uvetc=example3_setup(backward_in_time)
   else
      error("unknown example (nam parameter value)")
   end

   nx,ny=size(uvetc["XC"][1])

   lo0,lo1=lon_rng
   la0,la1=lat_rng
   #la0,la1=(45.0,55.0)
   n=1000
   lon=lo0 .+(lo1-lo0).*rand(n)
   lat=la0 .+(la1-la0).*rand(n)
   (u0,du)=initialize_lonlat(uvetc,lon,lat)
   u0[3,:] .= z_init

   #duvw(du,u0,uvetc,0.0)
   𝑇 = (0.0,uvetc["t1"]-uvetc["t0"])
   prob = ODEProblem(duvw,u0,𝑇,uvetc)
   #sol = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
   sol = solve(prob,Euler(),dt=uvetc["dt"])

   sol[1,:,:]=mod.(sol[1,:,:],nx)
   sol[2,:,:]=mod.(sol[2,:,:],ny)
   XC=exchange(uvetc["XC"])
   YC=exchange(uvetc["YC"])
   df=postprocess_lonlat(sol,XC,YC)

   #add third coordinate
   z=sol[3,:,:]
   df.z=z[:]

   #add time coordinate
   nt=size(df,1)/n
   df.t=uvetc["dt"]*[ceil(i/n)-1 for i in 1:nt*n]

   #to plot e.g. Pacific Ocean transports, shift longitude convention?
   df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0

   return df
end

"""
example3_setup()

Define gridded variables and return result as Dictionary (`uvetc`).
"""
function example3_setup(;backward_in_time::Bool=false)

   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/llc90_latlon/")

   γ=gcmgrid(dirIn,"PeriodicChannel",1,
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

   #mask out near edge values to avoid exiting domain
   u[:,1:2].=NaN
   v[:,1:2].=NaN
   u[:,end-2:end].=NaN
   v[:,end-2:end].=NaN

   u=read(u,MeshArray(γ,Float32))
   v=read(v,MeshArray(γ,Float32));

   u[findall(isnan.(u))]=0.0
   v[findall(isnan.(v))]=0.0

   backward_in_time ? s=-1.0 : s=1.0
   u0=s*u; u1=s*u;
   v0=s*v; v1=s*v;

   t0=0.0; t1=86400*366*2.0; dt=10*86400.0;

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


"""
OCCA_setup()

Define gridded variables and return result as Dictionary (`uvetc`).
"""
function OCCA_setup(;backward_in_time::Bool=false)

   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/GRID_LL360/")
   γ=GridSpec("PeriodicChannel",dirIn)
   Γ=GridLoad(γ)

   dirIn=joinpath(p,"../examples/OCCA_climatology/")
   n=length(Γ["RC"])

   fileIn=dirIn*"DDuvel.0406clim.nc"
   u = ncread(fileIn,"u")
   u=dropdims(mean(u,dims=4),dims=4)
   u[findall(u .< -1.0e10)] .=0.0
   u=read(u,MeshArray(γ,Float32,n))

   fileIn=dirIn*"DDvvel.0406clim.nc"
   v = ncread(fileIn,"v")
   v=dropdims(mean(v,dims=4),dims=4)
   v[findall(v .< -1.0e10)] .=0.0
   v=read(v,MeshArray(γ,Float32,n))

   fileIn=dirIn*"DDwvel.0406clim.nc"
   w = ncread(fileIn,"w")
   w=dropdims(mean(w,dims=4),dims=4)
   w[findall(w .< -1.0e10)] .=0.0
   w=-cat(w,zeros(360, 160),dims=3)
   w[:,:,1] .=0.0
   w=read(w,MeshArray(γ,Float32,n+1))

   for i in eachindex(u)
      u[i]=u[i]./Γ["DXC"][1]
      v[i]=v[i]./Γ["DYC"][1]
   end

   for i in eachindex(u)
      u[i]=circshift(u[i],[-180 0])
      v[i]=circshift(v[i],[-180 0])
   end

   for i in eachindex(w)
      w[i]=w[i]./Γ["DRC"][min(i[2]+1,n)]
      w[i]=circshift(w[i],[-180 0])
   end

   tmpx=circshift(Γ["XC"][1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ["XC"][1]=tmpx
   delete!.(Ref(Γ), ["hFacC", "hFacW", "hFacS"]);

   backward_in_time ? s=-1.0 : s=1.0
   u0=s*u; u1=s*u;
   v0=s*v; v1=s*v;
   w0=s*w; w1=s*w;

   t0=0.0; t1=86400*366*2.0; dt=dt=10*86400.0;

   uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1,
              "w0" => w0, "w1" => w1, "t0" => t0, "t1" => t1, "dt" => dt) ;

   return merge(uvt,Γ)

end
