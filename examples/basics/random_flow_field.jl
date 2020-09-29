# # Random Flow
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/random_flow_field.ipynb)
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
include(joinpath(p,"../examples/recipes_plots.jl"));

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

xy=transpose([x0[:] y0[:] ones(size(x0[:]))]);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Compute Trajectories

tr = DataFrame([fill(Int, 1) ; fill(Float64, 3)], [:ID, :x, :y, :t])
solv(prob) = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
𝐼 = Individuals{Float64}(📌=xy[:,:], 🔴=tr, 🆔=collect(1:size(xy,2)),
                         🚄 = dxy_dt!, ∫ = solv, 🔧 = postprocess_xy, 𝑃=𝑃)

𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.2 Plot Results
#
# For example, generate a simple animation (with `if true`):

if false
anim = @animate for t in 0:2.0:maximum(𝐼.tr.t)
   phi_and_subset(Γ,ϕ,𝐼.tr,t)
end
pth=tempdir()*"/"
gif(anim, pth*"RandomFlow.gif", fps = 15)
end
