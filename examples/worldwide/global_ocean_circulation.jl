
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
include(joinpath(p,"../examples/recipes_plots.jl"))
include(joinpath(p,"../examples/helper_functions.jl"))
get_grid_if_needed()

#velocity files:
pp="$p/../examples/nctiles_climatology"
q=dirname(pathof(OceanStateEstimation))
qq="$q/../examples/nctiles_climatology"
!isfile(pp*".csv") ? run(`cp $qq.csv $pp.csv`) : nothing
!isdir(pp) ? run(`mkdir $pp`) : nothing
!isdir(pp*"/UVELMASS") ? get_from_dataverse("UVELMASS",pp) : nothing
!isdir(pp*"/VVELMASS") ? get_from_dataverse("VVELMASS",pp) : nothing

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Set Up Parameters & Inputs
#
# - depth level and duration
# - read grid variables
# - read & normalize velocities

k=10 #choice of vertical level
ny=2 #number of simulated years (20 for k>20)
r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

#read grid and set up connections between subdomains
γ=GridSpec("LatLonCap",joinpath(p,"../examples/GRID_LLC90/"))
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

#initialize u0,u1 etc
𝑃=set_up_𝑃(k,0.0,Γ,joinpath(p,"../examples/nctiles_climatology/"));

#add parameters for use in reset!
tmp=(frac=r_reset, Γ=Γ)
𝑃=merge(𝑃,tmp)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# - initial particle positions randomly over Global Ocean
# - initial integration from time 0 to 0.5 month
# - update velocity fields & repeat for ny years

# ## 3.1 Initialize Individuals & Solution
#

println("Initialization started ...")

if false
    lo0, lo1 = (-160.0, -150.0)
    la0, la1 = (35.0, 45.0)
    np = 100
    lon = lo0 .+ (lo1 - lo0) .* rand(np)
    lat = la0 .+ (la1 - la0) .* rand(np)
    (u0, _) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
    id=collect(1:np)
else
    np=10000
    (lon, lat) = randn_lonlat(2*np)
    (u0, _) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
    u0=u0[:,1:np]
end

𝐼 = Individuals{Float64}(xy=deepcopy(u0), id=collect(1:np), 𝑃=deepcopy(𝑃))

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

[step!(𝐼) for y=1:ny, m=1:12]

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.4 Compute summary statistics
#

df = 𝐼.tr
gdf = groupby(𝐼.tr, :ID)
tmp = combine(gdf,nrow,:lat => mean)
show(tmp)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot trajectories
#

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Examples Using Various Plotting Packages are provided below.

#p=dirname(pathof(IndividualDisplacements))
#nn=1000

#include(joinpath(p,"../examples/recipes_plots.jl"))
#plt=PlotBasic(df,nn,180.)

#include(joinpath(p,"../examples/recipes_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

#include(joinpath(p,"../examples/recipes_makie.jl"))
#AbstractPlotting.inline!(true) #for Juno, set to false
#scene=PlotMakie(df,nn,180.0)
##Makie.save("LatLonCap300mDepth.png", scene)

#nb # %% {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# Here we create `lon`, `lat`, and `DL` (log10 of bottom depth) to use in plot background:

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))

DL=log10.(Interpolate(Γ["Depth"],f,i,j,w))
DL[findall((!isfinite).(DL))].=NaN
DL=reshape(DL,size(lon));

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Then we generate plot or movie using `GeoMakie.jl` ...

if false
#using ArgoData
#p = include(joinpath(dirname(pathof(ArgoData)),"movies.jl"));
df.y=2000 .+ df.t ./86400/365
yy=collect(2000:0.05:2000+ny)
scene = ProjMap(DL,colorrange=(2.,4.))
ProjScatterMovie(scene,df,yy,"GlobalDomain_fleet_k"*"$k"*"_v1.mp4",dt=1.0,mrksz=5e3)
end

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ... or using `Plots.jl`:

if false
contourf(lon[:,1],lat[1,:],transpose(DL),clims=(1.5,5),c = :ice, colorbar=false)

dt=0.0001
t0=0.0
t1=𝑃.𝑇[2]

t=t1
df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
scatter!(df_t.lon,df_t.lat,markersize=1.5,c=:red,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:circle, stroke(0)))

t=t0
df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
scatter!(df_t.lon,df_t.lat,markersize=1.5,c=:yellow,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:dot, stroke(0)))
end
