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

# + {"cell_style": "split", "cell_type": "markdown"}
# ## 1. Import Software
# -

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/plot_pyplot.jl"))

# ## 2. Setup Problem

# +
uvetc=IndividualDisplacements.example3_setup()

ii1=0:5:360; ii2=20:2:150;
n1=length(ii1); n2=length(ii2);
u0=Array{Float64,2}(undef,(2,n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2);
        i=i1+(i2-1)*n1
        u0[1,i]=ii1[i1]
        u0[2,i]=ii2[i2]       
end; end;

du=fill(0.0,size(u0))
⬡(du,u0,uvetc,0.0)
du
# -

# ## 3. Compute Trajectories
#
# - Define an ODE problem.
# - Solve the ODE problem to compute trajectories.

𝑇 = (0.0,uvetc["t1"]-uvetc["t0"])
prob = ODEProblem(⬡,u0,𝑇,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
size(sol)

# ## 4. Display results

# +
sol[1,:,:]=mod.(sol[1,:,:],360)
sol[2,:,:]=mod.(sol[2,:,:],180)
XC=exchange(uvetc["XC"])
YC=exchange(uvetc["YC"])
df=postprocess_ODESolution(sol,XC,YC)

PyPlot.figure(); PlotMapProj(df,3000)
# -


