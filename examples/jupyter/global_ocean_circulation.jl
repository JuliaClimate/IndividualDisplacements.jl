
# # Global Climatology
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/global_ocean_circulation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/global_ocean_circulation.ipynb)
#
# Advect particles with climatological monthly mean flow at selected depth level
# (e.g. `k=10` for 95 m) from a global ocean state estimate ([ECCO v4 r2](https://eccov4.readthedocs.io/en/latest/) ; see also <https://ecco-group.org>)
# which is here repeated for `ny` years. For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# [![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1. Get Software & Iput Files
#
# - packages + helper functions
# - grid and velocity files

using IndividualDisplacements, Statistics
import IndividualDisplacements.DataFrames: groupby, combine, nrow

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/worldwide/ECCO_FlowFields.jl"))

ECCO_FlowFields.get_ecco_velocity_if_needed()

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Set Up Parameters & Inputs
#
# - select vertical level (k=1 by default; k=0 for 3D) & duration in years (ny=2 by default)
# - read grid variables
# - return FlowFields (𝑃) and ancillary variables etc (𝐷) 
# - read & normalize velocities (𝐷.🔄)

𝑃,𝐷=ECCO_FlowFields.global_ocean_circulation(k=1,ny=2);

fieldnames(typeof(𝑃))

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# ### 3.1 Initialize Individuals & Solution
#
# - initial particle positions randomly over Global Ocean

np=10

#df = ECCO_FlowFields.init_global_randn(np,𝐷)
df = ECCO_FlowFields.init_from_file(np)

𝐼=Individuals(𝑃,df.x,df.y,df.f)
fieldnames(typeof(𝐼))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# - initial integration from time 0 to 0.5 month

📌ini=deepcopy(𝐼.📌)
𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### 3.2 Iteration function example
#
# In addition, `step!` is defined to provide additional flexibility around `∫!` :
#
# - `𝐷.🔄(𝐼.𝑃,t_ϵ)` resets the velocity input streams to bracket t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2]) 
# - `reset_📌!(𝐼)` randomly selects a fraction (`𝐷.frac`) of the particles and resets their positions before each integration period. This tends to maintain homogeneous coverage of the Global Ocean by particles.
# - `∫!(𝐼)` then solves for the individual trajectories over one month, with updated velocity fields (𝐼.𝑃.u0 etc), and adds diagnostics to the DataFrame used to record variables along the trajectory (𝐼.🔴).

function step!(𝐼::Individuals)
    t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2])
    𝐷.🔄(𝐼.𝑃,𝐷,t_ϵ)
    ECCO_FlowFields.reset_📌!(𝐼,𝐷.frac,📌ini)
    ∫!(𝐼)
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.3 Iterate For `ny*12` Months
#

[step!(𝐼) for y=1:1, m=1:1]

add_lonlat!(𝐼.🔴,𝐷.XC,𝐷.YC)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.4 Compute summary statistics
#
# See [DataFrames.jl](https://juliadata.github.io/DataFrames.jl/latest/) documentation for detail and additinal functionalities.

gdf = groupby(𝐼.🔴, :ID)
sgdf= combine(gdf,nrow,:lat => mean)
sgdf[rand(1:size(sgdf,1),4),:]

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot trajectories / individual positions
#
# ```
# using Plots
# p=plot(;xlims=(-180,180),ylims=(-90,90),legend=:none)
# p!(x,y)=scatter!(p,x,y,markersize=1.1,markerstrokewidth=0)
# [p!(gdf[i].lon,gdf[i].lat) for i in rand(collect(1:length(gdf)),10)]
# display(p)
# ```

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Or select a background map (e.g. `lon`, `lat`, and `DL=log10(bottom depth)`)
# and a recipe to superimpose initial and final locations. Try:
#
# ```
# include(joinpath(dirname(pathof(IndividualDisplacements)),"../examples/recipes_plots.jl"))
# map(𝐼,OceanDepthLog(𝐷.Γ))
# ```
