# # Three Dimensions
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/three_dimensional_ocean.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/three_dimensional_ocean.ipynb)
#
# Advect particles with climatological mean flow in three dimensions starting from a selected depth level
# (e.g. `k=10` for 95 m) and region using a near-global ocean state estimate ([OCCA](https://doi.org/10.1175/2009JPO4043.1)
# which is here repeated for two years. For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1. Load Software
#

using IndividualDisplacements, DataFrames, OceanStateEstimation
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/example123.jl"))
include(joinpath(p,"../examples/helper_functions.jl"))
get_ll360_grid_if_needed(); get_occa_velocity_if_needed();

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.1 Ocean Circulation Setup
#

nam="OCCA"
bck=false

if nam=="OCCA"
   𝑃,Γ=OCCA_setup(backward_in_time=bck)
   🚄 =dxyz_dt
elseif nam=="LL90"
   𝑃,Γ=example3_setup(backward_in_time=bck)
   🚄 =dxy_dt
else
   error("unknown example (nam parameter value)")
end

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.2 Solver And Analysis Setup
#

function ∫(prob)
   sol=solve(prob,Tsit5(),saveat=10*86400.0)
   nx,ny=𝑃.ioSize[1:2]
   sol[1,:,:]=mod.(sol[1,:,:],nx)
   sol[2,:,:]=mod.(sol[2,:,:],ny)
   return sol
end

function 🔧(sol,𝑃::NamedTuple;id=missing,𝑇=missing)
   df=postprocess_lonlat(sol,𝑃,id=id,𝑇=𝑇)

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=sol[3,:,:]
   df.k=k[:] #level
   k=Int.(floor.(df.k)); w=(df.k-k); 
   df.z=𝑃.RF[1 .+ k].*(1 .- w)+𝑃.RF[2 .+ k].*w #depth

   #add one isotherm depth
   θ=0.5*(𝑃.θ0+𝑃.θ1)
   d=isosurface(θ,15,𝑃.RC)
   d[findall(isnan.(d))].=0.
   df.iso=interp_to_xy(df,exchange(d));

   #add color = f(iso-z)
   c=fill(:gold,length(df.iso))
   c[findall(df.iso.<df.z)].=:violet
   df.col=c

   #to plot e.g. Pacific Ocean transports, shift longitude convention?
   df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
   return df
end

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.3 Initialize Individuals
#

"""
    set_up_individuals(𝑃,Γ,∫,🚄,🔧; nf=10000, z_init=4.5, 
               lon_rng=(-160.0,-150.0), lat_rng=(30.0,40.0))

Set up `Individuals` data structure with `nf` particles moving within a near-global Ocean domain. 
"""
function set_up_individuals(𝑃,Γ,∫,🚄,🔧; nf=10000, 
      z_init=4.5, lon_rng=(-160.0,-150.0), lat_rng=(30.0,40.0))

   lo0,lo1=lon_rng
   la0,la1=lat_rng

   lon=lo0 .+(lo1-lo0).*rand(nf)
   lat=la0 .+(la1-la0).*rand(nf)
   (xy,_)=initialize_lonlat(Γ,lon,lat)
   xy[3,:] .= z_init
   id=collect(1:size(xy,2))

   tr = DataFrame([fill(Int, 2) ; fill(Float64, 9); fill(Symbol, 1)], 
   [:ID, :fid, :x, :y, :k, :z, :iso, :t, :lon, :lat, :year, :col])

   𝐼 = Individuals{Float64}(📌=xy, 🔴=tr, 🆔=id, 🚄 = 🚄, ∫ = ∫, 🔧 = 🔧, 𝑃=𝑃)

   return 𝐼
end

set_up_individuals(𝐼::Individuals; nf=10000) = set_up_individuals(𝑃,Γ,∫,🚄,🔧; nf=nf)

𝐼=set_up_individuals(𝑃,Γ,∫,🚄,🔧,nf=100)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.1 Compute Displacements
#

𝑇=(0.0,𝐼.𝑃.𝑇[2])

∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.2 Analyze Results
#
# The recorded simulation output, 🔴, is a in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, which is easily manipulated or plotted.
#
# - either `Plots.jl`:

#md include(joinpath(p,"../examples/recipes_plots.jl"))
#md PlotBasic(𝐼.🔴,100,90.0)
#plot_end_points(𝐼,Γ)

# - or `Makie.jl`:

#include(joinpath(p,"../examples/recipes_Makie.jl"))
#p=PlotMakie(𝐼.🔴,100,180.);
#display(p)
