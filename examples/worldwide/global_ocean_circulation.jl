
# # Global Climatology
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/global_ocean_circulation.ipynb)
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

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Statistics, DataFrames, MITgcmTools, OceanStateEstimation

include(joinpath(dirname(pathof(IndividualDisplacements)),"../examples/helper_functions.jl"))
get_llc90_grid_if_needed(); get_ecco_velocity_if_needed();

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Set Up Parameters & Inputs
#
# - select vertical level & duration in years
# - read grid variables & velocities
# - normalize velocities

𝑃=setup_global_ocean(k=1,ny=2);

keys(𝑃)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# ### 3.1 Initialize Individuals & Solution
#
# - initial particle positions randomly over Global Ocean

xy=init_global_randn(1000,𝑃); id=collect(1:size(xy,2))
𝐼 = Individuals{Float64}(📌=xy[:,:], 🆔=id, 🔧=postprocess_lonlat, 🚄 = dxy_dt!, 𝑃=𝑃)

fieldnames(typeof(𝐼))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# - initial integration from time 0 to 0.5 month

𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### 3.2 Iteration function example
#
# - `𝐼.𝑃.🔄(𝐼.𝑃,t_ϵ)` resets the velocity input streams to bracket t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2]) 
# - `reset_lonlat!(𝐼)` randomly selects a fraction (defined in `setup_global_ocean()`) of the particles and resets their positions before each integration period. This can maintain homogeneous coverage of the Global Ocean by particles.
# - `∫!(𝐼,𝑇)` then solves for the individual trajectories over one month, after updating velocity fields (𝐼.u0 etc) if needed, and adds diagnostics to the DataFrame used to record / trace variables along the trajectory (𝐼.tr).

function step!(𝐼::Individuals)
    t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2])
    𝐼.𝑃.🔄(𝐼.𝑃,t_ϵ)
    reset_lonlat!(𝐼)
    𝑇=Tuple(𝐼.𝑃.𝑇)
    ∫!(𝐼,𝑇)
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.3 Iterate For `ny*12` Months
#

[step!(𝐼) for y=1:2, m=1:12]

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
# plot_end_points(𝐼,𝐼.𝑃.)
# ```
