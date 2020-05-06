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
using Plots, Statistics, MITgcmTools, DataFrames

# ## 2. Read gridded variables as `MeshArray`s

# Put grid variables in a dictionary.

γ=GridSpec("LatLonCap","GRID_LLC90/")
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))
uvetc=read_uvetc(20,Γ,"nctiles_climatology/");

# ## 3. Compute trajectories from gridded flow fields

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

plt=plot(u)
plot!(v)
display(plt)
# -

# Solve for trajectory in small test case.

𝑇 = (0.0,uvetc["t1"]-uvetc["t0"])
prob = ODEProblem(⬡!,uInit,𝑇,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
sol_two = solve(prob,Euler(),dt=1e6)
size(sol_one)

# Define initial condition array

(u0,du)=initialize_locations(uvetc,10);

# Solve for all trajectories.

prob = ODEProblem(⬡!,u0,𝑇,uvetc)
sol = solve(prob,Euler(),dt=uvetc["dt"])
size(sol)

# ## 4. Plot trajectories
#

# - Map i,j position to lon,lat coordinates and convert to DataFrame.

df=postprocess_ODESolution(sol,uvetc)

# - call `PlotMapProj`

# +
p=dirname(pathof(IndividualDisplacements))

nn=1000

include(joinpath(p,"../examples/plot_plots.jl"))
plt=PlotBasic(df,nn,180.)
display(plt)

#include(joinpath(p,"../examples/plot_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

#include(joinpath(p,"../examples/plot_makie.jl"))
#AbstractPlotting.inline!(true) #for Juno, set to false
#scene=PlotMakie(df,nn,180.0)
##Makie.save("LatLonCap300mDepth.png", scene)
# -
