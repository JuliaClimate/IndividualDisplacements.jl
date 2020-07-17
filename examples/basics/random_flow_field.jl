# # Random Flow Simulation
#
# Simulate trajectories of a particle cloud in a randomly generated flow field.
# A doubly periodic domain is used and an animation generated.
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/random_flow_field.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/random_flow_field.ipynb)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software

using OrdinaryDiffEq, IndividualDisplacements, MeshArrays
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))
p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"../examples/helper_functions.jl"))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.2 Setup Problem

# Put grid variables in a dictionary.

np=12
nq=24
Γ=simple_periodic_domain(np,nq);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Derive flow field from randomly generated ϕ streamfunction

𝑃,ϕ=setup_random_flow(Γ);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.3 Initial Conditions

x0=np*(0.1:0.02:0.9)
y0=nq*(0.1:0.02:0.9)

x0=vec(x0)*ones(1,length(y0))
y0=ones(size(x0,1),1)*transpose(vec(y0))

u0=transpose([x0[:] y0[:] ones(size(x0[:]))]);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Compute Trajectories

𝑇 = (𝑃["t0"],𝑃["t1"])
prob = ODEProblem(⬡!,u0,𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
size(sol)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.2 Process Output

df=postprocess_xy(sol,𝑃);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.3 Plot Results
#
# For example, to generate a simple animation:

include(joinpath(p,"../examples/recipes_plots.jl"))
anim = @animate for t in 0:2.0:maximum(df[!,:t])
   phi_and_subset(Γ,ϕ,df,t)
end
pth=tempdir()*"/"
gif(anim, pth*"RandomFlow.gif", fps = 15)
