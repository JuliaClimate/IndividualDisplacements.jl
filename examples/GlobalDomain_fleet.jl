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

# # This notebook
#
# _Notes:_ For documentation see <https://gaelforget.github.io/MeshArrays.jl/stable/>, <https://docs.juliadiffeq.org/latest/solvers/ode_solve.html> and <https://en.wikipedia.org/wiki/Displacement_(vector)>

# ## 1. import software

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Statistics, MITgcmTools, DataFrames

# ## 2. Read gridded variables as `MeshArray`s

# Put grid variables in a dictionary.

# +
k=1 #choice of vertical level
ny=10 #number of simulated years

#read grid and set up connections between subdomains
γ=GridSpec("LatLonCap","GRID_LLC90/")
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

#initialize u0,u1 etc
uvetc=IndividualDisplacements.read_uvetc(k,0.0,Γ,"nctiles_climatology/");
# -

# ## 3. Illustrate sample velocity & trajectory computations

# Let's illustrate the velocity interpolation scheme with a simple test first.

# +
uInit=[45.0,100.0,1.0]
du=fill(0.0,3);

ii=uInit[1]-3:0.1:uInit[1]+3
jj=uInit[2]-3:0.1:uInit[2]+3
ff=ones(size(jj))

s=size(ii)

(u,v,f)=[zeros(s),zeros(s),zeros(s)]
for i in eachindex(ii)
    ⬡!(du,[ii[i];jj[i];ff[i]],uvetc,0.0)
    u[i],v[i],f[i]=du
end

#using plots
#plt=plot(u)
#plot!(v)
#display(plt)
# -

# Solve for trajectory in small test case.

𝑇 = (0.0,uvetc["t1"])
prob = ODEProblem(⬡!,uInit,𝑇,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
sol_two = solve(prob,Euler(),dt=1e6)
size(sol_one)

# ## 3. Main computation of trajectories from gridded flow fields

# Define initial condition array

#(u0,du)=initialize_grid_locations(uvetc,10);
(u0,du)=initialize_random_locations(Γ,10000; msk=Γ["hFacC"][:,1]);

# Solve for all trajectories.

prob = ODEProblem(⬡!,u0,𝑇,uvetc)
sol = solve(prob,Euler(),dt=uvetc["dt"]/4.0)
size(sol)

# ## 4. Plot trajectories
#

# - Map i,j position to lon,lat coordinates and convert to DataFrame.

df=postprocess_ODESolution(sol,uvetc)
df[1:4,:]

t0=[uvetc["t1"]]
u0 = deepcopy(sol[:,:,end])
println(size(df))
for y=1:ny
    for m=1:12
        uvetc=IndividualDisplacements.read_uvetc(k,t0[1],Γ,"nctiles_climatology/");
        𝑇 = (uvetc["t0"],uvetc["t1"])
        prob = ODEProblem(⬡!,u0,𝑇,uvetc)
        sol = solve(prob,Euler(),dt=uvetc["dt"]/4.0)
        tmp = postprocess_ODESolution(sol[:,:,2:end],uvetc)
        append!(df,tmp)
        t0[1]=uvetc["t1"]
        u0[:,:] = deepcopy(sol[:,:,end])
    end
    println(size(df))
end

# - call `PlotBasic`, `PlotMapProj`, or `PlotMakie`

# +
#p=dirname(pathof(IndividualDisplacements))
#nn=1000

#include(joinpath(p,"../examples/plot_plots.jl"))
#plt=PlotBasic(df,nn,180.)
#display(plt)

#include(joinpath(p,"../examples/plot_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

#include(joinpath(p,"../examples/plot_makie.jl"))
#AbstractPlotting.inline!(true) #for Juno, set to false
#scene=PlotMakie(df,nn,180.0)
##Makie.save("LatLonCap300mDepth.png", scene)
# -

# - create `lon`, `lat`, and `DL` to use in plot background

# +
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
# -

# - generate movie using `GeoMakie.jl`

if true
    using ArgoData
    p = include(joinpath(dirname(pathof(ArgoData)),"movies.jl"));
    tt=collect(2000:0.05:2010)
    scene = ProjMap(DL,colorrange=(3.,4.))
    ProjScatterMovie(scene,df,tt,"tmp.mp4")
end
