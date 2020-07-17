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
# - <https://en.wikipedia.org/wiki/Displacement_(vector)>
# - <https://juliaclimate.github.io/IndividualDisplacements.jl/dev>
# - <https://juliaclimate.github.io/MeshArrays.jl/dev>
# - <https://docs.juliadiffeq.org/latest>

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

np=12
nq=24
Γ=simple_periodic_domain(np,nq);

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Derive flow field from randomly generated ϕ streamfunction
# -

𝑃,ϕ=setup_random_flow(Γ);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.3 Initial Conditions

# +
x0=np*(0.1:0.02:0.9)
y0=nq*(0.1:0.02:0.9)

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
# ## 2.2 Process Output
# -

df=postprocess_xy(sol,𝑃);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.3 Plot Results
#
# For example, to generate a simple animation:
#
# ```
# include("recipes_plots.jl")
# anim = @animate for t in 0:2.0:maximum(df[!,:t])
#    phi_and_subset(Γ,ϕ,df,t)
# end
# pth=tempdir()*"/"
# gif(anim, pth*"RandomFlow.gif", fps = 15)
# ```
