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

# ## 1. Import Software

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq, Statistics, DataFrames
p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"../examples/plot_plots.jl"))
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))

# ## 2. Setup Problem

# Put grid variables in a dictionary.

# +
γ,Γ=GridOfOnes("PeriodicDomain",16,20)
(_,ϕ,_,_)=demo2(Γ)

plt=heatmap(ϕ[1])
display(plt)
# -

# Derive random flow field using ϕ as a streamfunction

(u,v)=gradient(ϕ,Γ)
u=u./Γ["DXC"]#normalization to grid units
v=v./Γ["DYC"]
(u,v)=exchange(u,v,1)
u0=-v; u1=-v; v0=u; v1=u;

# Put velocity fields, time parameters, etc in a dictionary.

# +
uvt = Dict( "u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, 
            "t0" => 0.0, "t1" => 200.0, "dt" => 0.1)

uvetc=merge(uvt,Γ)#add grid variables

msk=Dict("mskW" => fill(1.0,u), "mskS" => fill(1.0,v))
uvetc=merge(uvetc,msk)

plt=heatmap(uvetc["mskW"][1,1].*uvetc["u0"][1,1],title="U at the start")
display(plt)

# +
ii1=0.5:0.5:20; ii2=0.5:0.5:20; ii3=[2.0 4.0 5.0 7.0 10.0 12.0 13.0 15.0];

n1=length(ii1); n2=length(ii2); n3=length(ii3);
u0=Array{Float64,2}(undef,(3,n3*n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2); for i3 in eachindex(ii3);
        i=i1+(i2-1)*n1+(i3-1)*n1*n2
        u0[1,i]=ii1[i1]
        u0[2,i]=ii2[i2]
        u0[3,i]=ii3[i3]
end; end; end;
# -

# ## 3. Compute Trajectories

𝑇 = (0.0,uvetc["t1"]-uvetc["t0"])
prob = ODEProblem(⬡!,u0,𝑇,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
size(sol)

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=sol[1,:,:]+20.0*mod.(sol[3,:,:].-1,4)
lat=sol[2,:,:]+20.0*floor.((sol[3,:,:].-1)/4)
fIndex=sol[3,:,:]
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:], fIndex=fIndex[:])
size(df)

nn=minimum([5000 size(u0,2)])
plt=PlotBasic(df,nn,10.0)
#savefig("PeriodicDomainRandomFlow.png")
display(plt)


