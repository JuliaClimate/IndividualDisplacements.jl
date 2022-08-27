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

using IndividualDisplacements, OceanStateEstimation
import IndividualDisplacements.DataFrames: DataFrame

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/jupyter/helper_functions.jl"))

OceanStateEstimation.get_occa_velocity_if_needed();

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.1 Ocean Circulation Setup
#

𝑃,𝐷,Γ=OCCA_FlowFields(nmax=5);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.3 Initialize Individual Positions
#

"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see the 
global ocean example.
"""
function initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])
   return x,y
end

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.3 Diagnostic / Post-Processing Setup
#

custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[], 
   k=Float64[], z=Float64[], iso=Float64[], t=Float64[], 
   lon=Float64[], lat=Float64[], year=Float64[], col=Symbol[])

function custom🔧(sol,𝑃::𝐹_MeshArray3D,𝐷::NamedTuple;id=missing,𝑇=missing)
   df=postprocess_MeshArray(sol,𝑃,𝐷,id=id,𝑇=𝑇)
   add_lonlat!(df,𝐷.XC,𝐷.YC)

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=[[sol[i][3,1] for i in 1:size(sol,3)];[sol[i][3,end] for i in 1:size(sol,3)]]
  
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

# ## 2.4 Individuals Data Structure
#
# Set up `Individuals` data structure with `nf` particles moving in 3D 
# on a regular 1 degree resolution grid covering most of the Globe.

nf=100; lo=(-160.0,-150.0); la=(30.0,40.0); kk=2.5; 
df=DataFrame(:z => fill(kk,nf),:f => fill(1,nf))
(df.x,df.y)=initial_positions(Γ, nf, lo, la)

𝐼=Individuals(𝑃,df.x,df.y,df.z,df.f,(🔴=custom🔴,🔧=custom🔧, 𝐷=𝐷))

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
# #p=plot(𝐼)
# p=plot_paths(𝐼.🔴,100,180.);
# display(p)
# ```

# ## 3.3 Alternatives (optional / unit testing)

(x,y,z,f)=𝐼.📌[1]
𝐽=Individuals(𝐼.𝑃,x,y,z,f)
diff(𝐼)
gcdist(𝐼);