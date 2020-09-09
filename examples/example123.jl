using MeshArrays, Statistics, OrdinaryDiffEq, IndividualDisplacements

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
Plots.plot(𝐼.tr.x,𝐼.tr.y,linewidth=5,lc=:black, title="One Trajectory Example",
xaxis="x",yaxis="y",label="Julia Solution") # legend=false
pl=Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,lc=:red,label="MITgcm Solution")
```
"""
function example2()
   𝑃,Γ=example2_setup()
   (xy,df,ref,nSteps)=example2_xy(𝑃)

   𝑃.𝑇[:] = [0.0,nSteps*3600.0]
   solv(prob) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
   tr = DataFrame( ID=[], x=[], y=[], t = [])

   𝐼 = Individuals{Float64}(xy=xy[:,:], 𝑃=𝑃, ⎔! = ⬡, □ = solv, ▽ = postprocess_xy, tr = tr)
   start!(𝐼)

   return 𝐼, df,ref
end

"""
example2_xy()

Read MITgcm/pkg/flt output
"""
function example2_xy(𝑃)
p=dirname(pathof(IndividualDisplacements))
dirIn=joinpath(p,"../examples/flt_example/")
prec=Float32
df=read_flt(dirIn,prec)
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
ref=ref./𝑃.dx
xy=[tmp[1,:lon];tmp[1,:lat]]./𝑃.dx
return xy,df,ref,nSteps
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

   ## Visualize velocity fields

   mskW=read_mds(γ.path*"hFacW",MeshArray(γ,Float32,nr))
   mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   mskS=read_mds(γ.path*"hFacS",MeshArray(γ,Float32,nr))
   mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)
   Γ=merge(Γ,Dict("mskW" => mskW, "mskS" => mskS))

   𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, dx=Γ["dx"],
        𝑇=[t0,t1], XC=Γ["XC"], YC=Γ["YC"], ioSize=(80,42))
   return 𝑃,Γ
end

"""
    example3(nam::String="OCCA" ; bck::Bool=false, z_init=0.5,
       lon_rng=(-165.0,-145.0), lat_rng=(25.0,35.0))

Run particle trajectory simulation over near-global ocean domain (79.5°S to
79.5°N in OCCA, or 69.5°S to 56.2°N in the regular-grid part of LLC90) and
return the result as a `DataFrame` that can be manipulated or plotted later.

```
using IndividualDisplacements, MAT, NetCDF, DataFrames, Plots
p=dirname(pathof(IndividualDisplacements))

𝐼=example3("OCCA");
df=𝐼.tr

include(joinpath(p,"../examples/recipes_plots.jl"))
PlotBasic(df,100,90.0)

#include(joinpath(p,"../examples/recipes_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,3000); gcf()

#include(joinpath(p,"../examples/recipes_Makie.jl"))
#PlotMakie(df,3000,180.)

nf=maximum(df.ID)
nt=size(df,1)/nf
dt=maximum(df.t)/(nt-1)
@gif for t in 0:nt-1
   scatter_zcolor(df,t*dt,df.z,(0,10))
end
```
"""
function example3(nam::String="OCCA" ; bck::Bool=false, z_init=0.5,
   lon_rng=(-165.0,-145.0), lat_rng=(25.0,35.0))
   if nam=="OCCA"
      𝑃,Γ=OCCA_setup(backward_in_time=bck)
   elseif nam=="LLC90"
      𝑃,Γ=example3_setup(backward_in_time=bck)
   else
      error("unknown example (nam parameter value)")
   end

   nx,ny=𝑃.ioSize[1:2]
   lo0,lo1=lon_rng
   la0,la1=lat_rng

   nf=100
   lon=lo0 .+(lo1-lo0).*rand(nf)
   lat=la0 .+(la1-la0).*rand(nf)
   (xy,_)=initialize_lonlat(Γ,lon,lat)
   xy[3,:] .= z_init
   id=collect(1:size(xy,2))

   function solv(prob)
      sol=solve(prob,Euler(),dt=10*86400.0)
      sol[1,:,:]=mod.(sol[1,:,:],nx)
      sol[2,:,:]=mod.(sol[2,:,:],ny)
      return sol
   end

   tr = DataFrame( ID=[], x=[], y=[], t = [], lon=[], lat=[], z=[], fid=[])

   function postproc(sol,𝑃::NamedTuple,id=missing)
      df=postprocess_lonlat(sol,𝑃,id)
      #add third coordinate
      z=sol[3,:,:]
      df.z=z[:]
      #to plot e.g. Pacific Ocean transports, shift longitude convention?
      df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
      return df
   end

   𝐼 = Individuals{Float64}(xy=xy, id=id, 𝑃=𝑃, ⎔! = dxyz_dt, □ = solv, ▽ = postproc, tr = tr)
   start!(𝐼)

   return 𝐼
end

"""
    example3(lon_rng,lat_rng,z_init,backward_in_time)

Run particle trajectory simulation over near-global ocean domain (79.5°S to
79.5°N in OCCA) and return the result as a `gif` animation using `scatter_zcolor`
and `Plots.jl`

```
using IndividualDisplacements, NetCDF, Plots
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))

example3((-165.0,-155.0),(5.0,15.0),5.5,true)
```
"""
function example3(lon_rng,lat_rng,z_init,bck)
   df=example3("OCCA",bck=bck, z_init=z_init,lon_rng=lon_rng,lat_rng=lat_rng)
   nf=maximum(df.ID)
   nt=size(df,1)/nf
   dt=maximum(df.t)/(nt-1)
   return @gif for t in 0:nt-1
        scatter_zcolor(df,t*dt,df.z,(0,10))
   end
end

"""
    example3_setup(;backward_in_time::Bool=false)

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

   t0=0.0; t1=86400*366*2.0;

   u0=u0./Γ["DXC"]
   u1=u1./Γ["DXC"]
   v0=v0./Γ["DYC"]
   v1=v1./Γ["DYC"]

   nr=50; kk=1;

   mskW=read(γ.path*"hFacW.latlon.data",MeshArray(γ,Float32,nr))
   mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   mskS=read(γ.path*"hFacS.latlon.data",MeshArray(γ,Float32,nr))
   mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)
   msk=Dict("mskW" => mskW, "mskS" => mskS)

   𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, 𝑇=[t0,t1], ioSize=(360,178),
        XC=exchange(Γ["XC"]), YC=exchange(Γ["YC"]))

   return 𝑃,Γ

end

"""
    OCCA_setup(;backward_in_time::Bool=false)

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

   tmpx=circshift(Γ["XG"][1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ["XG"][1]=tmpx
   Γ["Depth"][1]=circshift(Γ["Depth"][1],[-180 0])

   delete!.(Ref(Γ), ["hFacC", "hFacW", "hFacS","DXG","DYG","RAC","RAZ","RAS"]);

   backward_in_time ? s=-1.0 : s=1.0
   u0=s*u; u1=s*u;
   v0=s*v; v1=s*v;
   w0=s*w; w1=s*w;

   t0=0.0; t1=86400*366*2.0;

   𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, w0=w0, w1=w1, 𝑇=[t0,t1],
   XC=exchange(Γ["XC"]), YC=exchange(Γ["YC"]), ioSize=(360,160,50))

   return 𝑃,Γ

end
