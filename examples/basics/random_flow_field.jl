# # Random Flow Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/random_flow_field.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/random_flow_field.ipynb)
#
# Simulate trajectories of a particle cloud in a randomly generated flow field.
# A doubly periodic domain is used and an animation generated.
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# ![particles in random flow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software

using OrdinaryDiffEq, IndividualDisplacements, MeshArrays
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
include(joinpath(p,"../examples/recipes_plots.jl"))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.2 Setup Problem

# Put grid variables in a dictionary.

np=8
nq=12
Γ=simple_periodic_domain(np,nq);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Derive flow field from randomly generated ϕ streamfunction

𝑃,ϕ=setup_random_flow(Γ);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.3 Initial Conditions

x0,x1=np .*(0.4,0.6)
y0,y1=np .*(0.4,0.6)

n=100
x0=x0 .+(x1-x0).*rand(n)
y0=y0 .+(y1-y0).*rand(n)

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
# For example, generate a simple animation (with `if true`):

if false
anim = @animate for t in 0:2.0:maximum(df[!,:t])
   phi_and_subset(Γ,ϕ,df,t)
end
pth=tempdir()*"/"
gif(anim, pth*"RandomFlow.gif", fps = 15)
end
