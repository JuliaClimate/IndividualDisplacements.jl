
# # Global Ocean Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/global_ocean_circulation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/global_ocean_circulation.ipynb)
#
# Particles moving with climatological monthly mean flow at selected depth level
# (e.g. `k=10` for 95 m) based on an ocean state estimate (ECCO v4 r2 from https://ecco-group.org).
# For additional documentation e.g. see :
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

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
get_grid_if_needed()
get_velocity_if_needed()

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Set Up Parameters & Inputs
#
# - depth level and duration
# - read grid variables
# - read & normalize velocities

𝑃=setup_global_ocean(10);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# - initial particle positions randomly over Global Ocean
# - initial integration from time 0 to 0.5 month
# - update velocity fields & repeat for ny years

# ## 3.1 Initialize Individuals & Solution
#

u0=init_global_randn(10000,𝑃)
np=size(u0,2)
𝐼 = Individuals{Float64}(xy=deepcopy(u0), id=collect(1:np), 𝑃=deepcopy(𝑃))
fieldnames(typeof(𝐼))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Solve for displacements over first 1/2 month

start!(𝐼)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.2 Define iteration function
#
# _A fraction of the particles, randomly selected, is reset every month to maintain a relatively homogeneous coverage of the Global Ocean by the fleet of particles._

#𝐼 = Individuals{Float64}(xy=deepcopy(u0), id=deepcopy(id), 𝑃=deepcopy(𝑃), tr=deepcopy(df))

function step!(𝐼::Individuals)
    reset!(𝐼)
    displace!(𝐼)
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.3 Iterate For `ny*12` Months
#

println("Main loop started ...")

[step!(𝐼) for y=1:2, m=1:12]

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.4 Compute summary statistics
#

gdf = groupby(𝐼.tr, :ID)
tmp = combine(gdf,nrow,:lat => mean)
show(tmp)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot trajectories / individual positions
#

include(joinpath(p,"../examples/recipes_plots.jl"))
nn=min(length(𝐼.id),100)
plt=PlotBasic(𝐼.tr,nn,180.)

#nb # %% {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# Here we create `lon`, `lat`, and `DL` (log10 of bottom depth) to use in plot background:

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,_,_,_)=InterpolationFactors(𝐼.𝑃.Γ,vec(lon),vec(lat))

DL=log10.(Interpolate(𝐼.𝑃.Γ["Depth"],f,i,j,w))
DL[findall((!isfinite).(DL))].=NaN
DL=reshape(DL,size(lon));

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Then we generate plot or movie using e.g. `Plots.jl`:

contourf(lon[:,1],lat[1,:],transpose(DL),clims=(1.5,5),c = :ice, colorbar=false)

dt=0.0001
t0=0.0
t1=𝑃.𝑇[2]

t=t1
df = 𝐼.tr[ (𝐼.tr.t.>t-dt).&(𝐼.tr.t.<=t) , :]
scatter!(df.lon,df.lat,markersize=1.5,c=:red,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:circle, stroke(0)))

t=t0
df = 𝐼.tr[ (𝐼.tr.t.>t-dt).&(𝐼.tr.t.<=t) , :]
scatter!(df.lon,df.lat,markersize=1.5,c=:yellow,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:dot, stroke(0)))
