# # Particle Set
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/random_flow_field.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/random_flow_field.ipynb)
#
# Simulate trajectories of a particle cloud in a two-dimensional flow field.
# A doubly-periodic domain and randomly-generated flow fields are initially used.
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# Exercises: 
# - change the initial distribution of partices
# - increase the duration of the trajectories simulation
# - treat the non-periodic domain case by padding `u,v` with zeros 
# - replace `u,v` with your own two-dimensional flow fields 
#
# ![particles in random flow](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/RandomFlow.gif)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1. Import Software

using IndividualDisplacements, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"));

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2. Flow Fields
#

u,v,ϕ=random_flow_field();

# The above `u,v` arrays can be replaced with any other pair provided by the user.
#
# A couple of important considerations, however:
#
# - `u,v` are staggered on a C-grid; by `-0.5` grid point in direction `1` for `u` (`2` for `v`)
#  from the grid cell center (0.5,0.5)
# - `u,v` here derive from streamfunction `ϕ`, defined at the corner point, which ensures that 
#  the resulting `u,v` is non-divergent, purely rotational, over the C-grid domain.
# In brief:
#
# ```
# u=-(circshift(ϕ, (0,-1))-ϕ)
# v=(circshift(ϕ, (-1,0))-ϕ)
# ```

# If user were to start with collocated velocity (`uC,vC` at the grid cell center) then
# one can easily obtain the staggered velocity (`u,v`) as follows. These may contain both 
# [rotational and divergent](https://en.wikipedia.org/wiki/Helmholtz_decomposition) components.
#
# ```
# u=0.5*(circshift(uC, (0,1))+uC)
# v=0.5*(circshift(vC, (1,0))+vC)
# ```

# A convenient way to set up the flow fields using the MeshArrays.jl package (which 
# handles such staggered grids in general fashion) is to call `convert_to_FlowFields`

𝐹=convert_to_FlowFields(u,v);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Initialize Individuals
#
# For example, we can initialize 100 particles within a central subdomain as follows.

np,nq=size(u)
x=np*(0. .+ 1.0*rand(1000))
y=nq*(0. .+ 1.0*rand(1000))
a=ones(size(x)); #subdomain array index (just 1 here)

# The following constructor function wraps everything in the `Individuals` data structure.

𝐼=Individuals(𝐹,x,y,a)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Compute Trajectories
#
# The time period is `𝐼.𝑃.𝑇` by default, unless `∫!(𝐼,𝑇)` is called instead.

∫!(𝐼)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot Results
#
# For example, generate a simple animation:

#md p=dirname(pathof(IndividualDisplacements))
#md include(joinpath(p,"../examples/recipes_plots.jl"));

#md 🔴_by_t = groupby(𝐼.🔴, :t)
#md anim = @animate for t in eachindex(🔴_by_t)
#md    phi_scatter(ϕ,🔴_by_t[t])
#md end

#md pth=tempdir()*"/"
#md gif(anim, pth*"RandomFlow.gif", fps = 10)
