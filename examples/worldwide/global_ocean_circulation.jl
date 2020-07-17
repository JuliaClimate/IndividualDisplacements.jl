# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# # Global Ocean Particle Simulation Example
#
# Particles moving with climatological monthly mean flow at selected depth level (e.g. `k=10` for 95 m) based on an ocean state estimate (ECCO v4 r2 from https://ecco-group.org).

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### For More Documentation 
#
# - <https://en.wikipedia.org/wiki/Displacement_(vector)>
# - <https://juliaclimate.github.io/IndividualDisplacements.jl/dev>
# - <https://juliaclimate.github.io/MeshArrays.jl/dev>
# - <https://docs.juliadiffeq.org/latest> 

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1. import software

# +
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Statistics, MITgcmTools, DataFrames

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
get_grid_if_needed()

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Read gridded variables as `MeshArray`s

# +
k=10 #choice of vertical level
ny=10 #number of simulated years (20 for k>20)
r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

#read grid and set up connections between subdomains
γ=GridSpec("LatLonCap","GRID_LLC90/")
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

#initialize u0,u1 etc
uvetc=read_uvetc(k,0.0,Γ,"nctiles_climatology/");

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3. Sample velocity & trajectory computations

# +
#initial condition
uInit=[45.0,100.0,1.0]
du=fill(0.0,3);

#velocity computation over neighbor sample
ii=uInit[1]-3:0.1:uInit[1]+3
jj=uInit[2]-3:0.1:uInit[2]+3
ff=ones(size(jj))

s=size(ii)
(u,v,f)=[zeros(s),zeros(s),zeros(s)]
for i in eachindex(ii)
    ⬡!(du,[ii[i];jj[i];ff[i]],uvetc,0.0)
    u[i],v[i],f[i]=du
end

# + {"slideshow": {"slide_type": "subslide"}}
#using plots
#plt=plot(u)
#plot!(v)
#display(plt)

# + {"slideshow": {"slide_type": "subslide"}}
#Solve for single trajectory (uInit)
𝑇 = (0.0,uvetc["t1"])
prob = ODEProblem(⬡!,uInit,𝑇,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
sol_two = solve(prob,Euler(),dt=1e6)
size(sol_one)

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# - initial particle positions randomly over Global Ocean
# - initial integration from time 0 to 0.5 month
# - update velocity fields & repeat for n years
# -

# ## 3.1 Initialization & Initial Solution  

# +
#(u0,du)=initialize_gridded(uvetc,10);

#(lon, lat) = randn_lonlat(20000)
#(u0,du)=initialize_lonlat(Γ,lon,lat; msk=Γ["hFacC"][:,k]);

lo0,lo1=(-160.0,-150.0)
la0,la1=(35.0,45.0)
n=1000
lon=lo0 .+(lo1-lo0).*rand(n)
lat=la0 .+(la1-la0).*rand(n)
(u0,du)=initialize_lonlat(Γ,lon,lat; msk=Γ["hFacC"][:,k]);

# + {"slideshow": {"slide_type": "skip"}}
u0_store = deepcopy(u0)
n_store = size(u0_store,2)

#r_reset = 0.05 #fraction of the particles reset per month
n_reset = Int(round(r_reset*n_store))
#k_reset = rand(1:size(u0_store,2), n_reset)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Solve for all trajectories for first 1/2 month
# -

prob = ODEProblem(⬡!,u0,𝑇,uvetc)
sol = solve(prob,Euler(),dt=uvetc["dt"]/8.0)
size(sol)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Map `i,j` to `lon,lat` coordinates and convert to `DataFrame`
# -

df=postprocess_lonlat(sol,uvetc)
df[1:4,:]

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.2 Repeat `ny` Years x 12 Months
#
# _A fraction of the particles, randomly selected, is reset every month to maintain a relatively homogeneous coverage of the Global Ocean by the fleet of particles._
# -

t0=[uvetc["t1"]]
u0 = deepcopy(sol[:,:,end])
println(size(df))
for y=1:ny
    for m=1:12
        uvetc=read_uvetc(k,t0[1],Γ,"nctiles_climatology/")
        𝑇 = (uvetc["t0"],uvetc["t1"])
        prob = ODEProblem(⬡!,u0,𝑇,uvetc)
        sol = solve(prob,Euler(),dt=uvetc["dt"]/8.0)
        tmp = postprocess_lonlat(sol[:,:,2:end],uvetc)

        k_reset = rand(1:size(u0_store,2), n_reset)
        k_new = rand(1:size(u0_store,2), n_reset)
        t_reset = Int(size(tmp,1)/n_store)-1

        tmp[k_reset.+t_reset*n_store,2:end].=NaN #reset a random subset of particles
        append!(df,tmp)

        t0[1]=uvetc["t1"]
        u0[:,:] = deepcopy(sol[:,:,end])
        u0[:,k_reset].=deepcopy(u0_store[:,k_new]) #reset a random subset of particles
    end
    println(size(df))
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot trajectories
#

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# - Try `PlotBasic`, `PlotMapProj`, or `PlotMakie`

# +
#p=dirname(pathof(IndividualDisplacements))
#nn=1000

#include(joinpath(p,"../examples/recipes_plots.jl"))
#plt=PlotBasic(df,nn,180.)
#display(plt)

#include(joinpath(p,"../examples/recipes_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

#include(joinpath(p,"../examples/recipes_makie.jl"))
#AbstractPlotting.inline!(true) #for Juno, set to false
#scene=PlotMakie(df,nn,180.0)
##Makie.save("LatLonCap300mDepth.png", scene)

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# - create `lon`, `lat`, and `DL` to use in plot background and ...

# + {"slideshow": {"slide_type": "skip"}}
nf=size(u0,2)
nt=size(df,1)/nf
t=[ceil(i/nf)-1 for i in 1:nt*nf]
df[!,:t]=2000 .+ uvetc["dt"]/4/86400/365 * t

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))

DL=log10.(Interpolate(Γ["Depth"],f,i,j,w))
DL[findall((!isfinite).(DL))].=NaN
DL=reshape(DL,size(lon));

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# - generate movie using `GeoMakie.jl`
# -

if false
    using ArgoData
    p = include(joinpath(dirname(pathof(ArgoData)),"movies.jl"));
    tt=collect(2000:0.05:2000+ny)
    scene = ProjMap(DL,colorrange=(2.,4.))
    ProjScatterMovie(scene,df,tt,"GlobalDomain_fleet_k"*"$k"*"_v1.mp4",dt=1.0,mrksz=5e3)
end

# +
if false
    using Plots
    dt=0.0001
    t0=minimum(df[!,:t])
    t1=maximum(df[!,:t])

    #t=2001.0
    t=t1
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    scatter(df_t.lon,df_t.lat,markersize=1.5,c=:red,leg=:none,
        xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:circle, stroke(0)))
    t=t0
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    scatter!(df_t.lon,df_t.lat,markersize=0.1,c=:blue,leg=:none,
        xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:dot, stroke(0)))

end
# -


