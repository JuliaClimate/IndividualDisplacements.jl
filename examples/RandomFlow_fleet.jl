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

# # Particle Cloud Simulation Example
#
# In this example we simulate the trajectory of a particle cloud in a randomly generated flow field, in a doubly periodic domain. As an exercise left to the user, directions are provided e.g. to modify the size of the domain or the rate of divergence within the particle cloud.
#
# - 1. setup the software and initialize example
# - 2. simulate trajectories & plot results
# - 3. experiment with parameters (user)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### For More Documentation
#
# - <https://docs.juliadiffeq.org/latest>
# - <https://en.wikipedia.org/wiki/Displacement_(vector)>
# - <https://juliaclimate.github.io/IndividualDisplacements.jl/dev>
# - <https://juliaclimate.github.io/MeshArrays.jl/dev>

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software
# -

using OrdinaryDiffEq, IndividualDisplacements, MeshArrays
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))
p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"../examples/helper_functions.jl"))

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.2 Setup Problem
# -

# Put grid variables in a dictionary.

# +
np=24

Γ=SetupPeriodicDomain(np);

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Derive flow field from randomly generated ϕ streamfunction
# -

𝑃,ϕ=SetupRandomFlow(Γ);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.3 Initial Conditions

# +
x0=np*(0.2:0.005:0.25)
y0=np*(0.6:0.005:0.65)

x0=vec(x0)*ones(1,length(y0))
y0=ones(size(x0,1),1)*transpose(vec(y0))

u0=transpose([x0[:] y0[:] ones(size(x0[:]))]);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Compute Trajectories
# -

𝑇 = (𝑃["t0"],𝑃["t1"])
prob = ODEProblem(⬡!,u0,𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
size(sol)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.2 Process Output For Plotting
# -

using DataFrames

#x,y axes etc
x=sol[1,:,:]
y=sol[2,:,:]
fIndex=sol[3,:,:]
ID=collect(1:size(sol,2))*ones(1,size(sol,3))
df = DataFrame(ID=Int.(ID[:]), x=mod.(x[:],Ref(np)), y=mod.(y[:],Ref(np)), fIndex=fIndex[:]);

#time axis
nf=size(u0,2)
nt=size(df,1)/nf
t=[ceil(i/nf)-1 for i in 1:nt*nf]
df[!,:t]=(𝑃["t1"]-𝑃["t0"])/t[end].*t;

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.3 Plot Results

# +
using Plots, ColorSchemes

function scatter_subset(Γ,df,t,dt=5.0)
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    contourf(vec(Γ["XC"][1][:,1]),vec(Γ["YC"][1][1,:]),transpose(ϕ[1]),c = :blues,linewidth = 0.1)
    scatter!(df_t.x,df_t.y,markersize=2.0,c=:red,
    xlims=(0,np),ylims=(0,np),leg=:none,marker = (:circle, stroke(0)))
end

anim = @animate for t in 0:2.0:maximum(df[!,:t])
   scatter_subset(Γ,df,t)
end
pth=tempdir()*"/"
gif(anim, pth*"RandomFlow.gif", fps = 15)
# -


