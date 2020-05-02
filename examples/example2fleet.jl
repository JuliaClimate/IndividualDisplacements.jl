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
#     display_name: Julia 1.3.0-rc4
#     language: julia
#     name: julia-1.3
# ---

# # This notebook
#
# _Notes:_ For documentation see <https://gaelforget.github.io/MeshArrays.jl/stable/>, <https://docs.juliadiffeq.org/latest/solvers/ode_solve.html> and <https://en.wikipedia.org/wiki/Displacement_(vector)>

# ## 1. import software

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Plots, Statistics, MITgcmTools, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"plot_Plots.jl"))

# ## 2. Read gridded variables as `MeshArray`s

dirIn="flt_example/"
prec=Float32
df=ReadDisplacements(dirIn,prec)
uvetc=IndividualDisplacements.example2_setup()

# ## 2. Recompute displacements from gridded flow fields

# Define the ODE problem

# +
nSteps=2998
tspan = (0.0,nSteps*3600.0)

comp_vel=IndividualDisplacements.VelComp
get_vel=IndividualDisplacements.VelCopy
# -

# Set up initial conditions
#
# _Note how the size of sol (ie nb of steps) depends on initial location:_

# +
#ii1=1:10:80; ii2=1:10:42; #->sol is (2, 40, 40065)
#ii1=30:37; ii2=16:20; #->sol is (2, 40, 9674)
#ii1=10:17; ii2=16:20; #->sol is (2, 40, 51709)
ii1=5:5:40; ii2=5:5:25; #->sol is (2, 40, 51709)

n1=length(ii1); n2=length(ii2);
uInitS=Array{Float64,2}(undef,(2,n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2);
        i=i1+(i2-1)*n1
        uInitS[1,i]=ii1[i1]-0.5
        uInitS[2,i]=ii2[i2]-0.5
end; end;

prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
# -

# Compute solution

sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
size(sol)

# Display results

# +
ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=5000* mod.(sol[1,:,:],80)
lat=5000* mod.(sol[2,:,:],42)
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:])

plt=PlotBasic(df,size(sol,2),100000.0)
display(plt)
# -
