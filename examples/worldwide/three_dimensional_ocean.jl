# # Three Dimensions
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/three_dimensional_ocean.ipynb)
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
# ![Three dimensional simulation](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1. Load Software
#

using IndividualDisplacements, DataFrames

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
IndividualDisplacements.get_occa_velocity_if_needed();

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.1 Ocean Circulation Setup
#

𝑃,𝐷,Γ=OCCA_FlowFields(backward_in_time=false);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.2 Post-Processor Setup
#

tr = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[], 
               k=Float64[], z=Float64[], iso=Float64[], t=Float64[], 
               lon=Float64[], lat=Float64[], year=Float64[], col=Symbol[])

function 🔧(sol,𝑃::𝐹_MeshArray3D;id=missing,𝑇=missing)
   df=postprocess_MeshArray(sol,𝑃,id=id,𝑇=𝑇)
   add_lonlat!(df,𝐷.XC,𝐷.YC)

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=[sol[1,i,j][3] for i in 1:size(sol,2), j in 1:size(sol,3)]
   nz=length(𝐼.𝑃.u1)
   df.k=min.(max.(k[:],Ref(0.0)),Ref(nz)) #level
   k=Int.(floor.(df.k)); w=(df.k-k); 
   df.z=𝐷.RF[1 .+ k].*(1 .- w)+𝐷.RF[2 .+ k].*w #depth

   #add one isotherm depth
   θ=0.5*(𝐷.θ0+𝐷.θ1)
   d=isosurface(θ,15,𝐷.RC)
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
# ## 2.3 Initialize Individual Positions
#

"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see 
`initialize_lonlat()`.
"""
function initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

   lo0,lo1=lon_rng
   la0,la1=lat_rng

   lon=lo0 .+(lo1-lo0).*rand(nf)
   lat=la0 .+(la1-la0).*rand(nf)
   x=lon .+ (21. - Γ["XC"][1][21,1])
   y=lat .+ (111. - Γ["YC"][1][1,111])

   return x,y
end

# ## 2.4 Individuals data structure
#

"""
    set_up_individuals(𝑃,Γ,🔧; nf=10000, z_init=4.5, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Set up `Individuals` data structure with `nf` particles moving 
on a regular 1 degree resolution grid covering most of the Globe.
"""

function set_up_individuals(𝑃,Γ,🔧; nf=10000,
   z_init=2.5, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

   (x,y)=initial_positions(Γ; nf, lon_rng, lat_rng)
   xy = permutedims([Float32.([x[i];y[i];z_init;1.0]) for i in eachindex(x)])
   I=(position=xy,record=deepcopy(tr),velocity=dxyz_dt!,postprocessing=🔧,parameters=𝑃)
   𝐼=Individuals(I)
   return 𝐼
end

set_up_individuals(𝐼::Individuals; nf=10000) = set_up_individuals(𝑃,Γ,🔧; nf=nf)

𝐼=set_up_individuals(𝑃,Γ,🔧,nf=100)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.1 Compute Displacements
#

𝑇=(0.0,10*86400.0)

∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.2 Analyze Results
#
# The recorded simulation output, 🔴, is a in the [DataFrames](https://juliadata.github.io/DataFrames.jl/latest/) tabular format, which is easily manipulated or plotted.

# - either `Plots.jl`:
#
# ```
# include(joinpath(p,"../examples/recipes_plots.jl"))
# p=plot(𝐼)
# #p=map(𝐼,OceanDepthLog(Γ))
# display(p)
# ```

# - or `Makie.jl`:
#
# ```
# include(joinpath(p,"../examples/recipes_Makie.jl"))
# p=PlotMakie(𝐼.🔴,100,180.);
# display(p)
# ```
