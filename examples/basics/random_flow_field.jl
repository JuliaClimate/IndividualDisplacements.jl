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
# ## 1. Import Software

using IndividualDisplacements, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Setup Problem

# ### 2.1 Sample flow field
#
# The `u,v` arrays below can be replaced with any other pair provided by the user.
#
# A couple of important considerations, however:
# - `u,v` are staggered on a C-grid; by `-0.5` grid point in direction `1` for `u` (`2` for `v`)
#  from the grid cell center (0.5,0.5)
# - `u,v` here derive from streamfunction `ϕ`, defined at the corner point, which ensures that 
#  the resulting `u,v` is non-divergent, purely rotational, over the C-grid domain
#
# In brief:
# ```
# u=-(circshift(ϕ, (0,-1))-ϕ)
# v=(circshift(ϕ, (-1,0))-ϕ)
# ```
#
# If user were to start with collocated velocity (`uC,vC` at the grid cell center) then
# one can easily obtain the staggered velocity (`u,v`) as follows. These may contain both 
# [rotational and divergent](https://en.wikipedia.org/wiki/Helmholtz_decomposition) components.
#
# ```
# u=0.5*(circshift(uC, (0,1))+uC)
# v=0.5*(circshift(vC, (1,0))+vC)
# ```

u,v,ϕ=setup_random_flow()

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 2.2 Initialize Individuals

np,nq=size(u)
x=np*(0.4 .+ 0.2*rand(100))
y=nq*(0.4 .+ 0.2*rand(100))

𝐼=setup_point_cloud(u,v,X=x,Y=y)
#𝐼.𝑃.𝑇[2]=1000.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Compute Trajectories

∫!(𝐼)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot Results
#
# For example, generate a simple animation:

#!jl p=dirname(pathof(IndividualDisplacements))
#!jl include(joinpath(p,"../examples/recipes_plots.jl"));

#!jl 🔴_by_t = groupby(𝐼.🔴, :t)
#!jl anim = @animate for t in eachindex(🔴_by_t)
#!jl    phi_scatter(ϕ,🔴_by_t[t])
#!jl end

#!jl pth=tempdir()*"/"
#!jl gif(anim, pth*"RandomFlow.gif", fps = 15)
