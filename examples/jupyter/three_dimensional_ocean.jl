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

using IndividualDisplacements
import IndividualDisplacements.DataFrames: DataFrame

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/worldwide/OCCA_FlowFields.jl"))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.1 Ocean Circulation Setup
#

𝑃,𝐷=OCCA_FlowFields.setup(nmax=5);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2.2 Initialize Positions
#

nf=100; lo=(-160.0,-150.0); la=(30.0,40.0); level=2.5; 
df=OCCA_FlowFields.initial_positions(𝐷.Γ, nf, lo, la, level)

# ## 2.4 Individuals Data Structure
#
# Set up `Individuals` data structure with `nf` iindividuals moving in 3D 
# on a regular 1 degree resolution grid covering most of the Globe.

𝐼=Individuals(𝑃,df.x,df.y,df.z,df.f,
   (🔴=OCCA_FlowFields.custom🔴,🔧=OCCA_FlowFields.custom🔧, 𝐷=𝐷))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 3.1 Compute Displacements
#

𝑇=(0.0,10*86400.0)

∫!(𝐼,𝑇)

# ## 3.2 Alternatives (optional / unit testing)

(x,y,z,f)=𝐼.📌[1]
𝐽=Individuals(𝐼.𝑃,x,y,z,f)
diff(𝐼)
gcdist(𝐼);