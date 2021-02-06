# # Three Dimensions
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/solid_body_rotation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/solid_body_rotation.ipynb)
#
# Simulate the trajectory of an individual point, first in a perfectly circular flow (a.k.a. solid body rotation). Then add a convergent term to obtain a spiraling trajectory, and a constant vertical velocity for the third dimension. These simple flow configurations can be thought of as idealized models e.g. ocean meso-scale eddies.
#
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# Exercise examples: 
# - make the sinking velocity decrease with time (hint: it increases in the original notebook) 
# - change the number of times the particle goes around the origin
#
# ![solid body rotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# # 1 Problem Configuration
#
# Here we set up software, grid, flow fields, initial conditions.
#
# ### 1.1 Import Software

using IndividualDisplacements, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/flow_fields.jl"))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### 1.2  Flow Fields
#
# The `simple_flow_field` function (defined in `helper_functions.jl`) defines a simple
# three-dimensional flow field. Exercise: locate `simple_flow_field` and modify the 
# flow field parameters (e.g. intensity and sign of the convergent term).

np,nz=16,4 #gridded domain size (horizontal and vertical)

u,v,w=solid_body_rotation(np,nz) #staggered velocity arrays

𝐹=𝐹_Array3D{eltype(u)}(u,u,v,v,0*w,1*w,[0,19.95*2*pi]); #FlowFields data structure

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 1.3 Initialize Individuals
#
# Let's just set up one individual at [np*1/3,np*1/3,nz*1/3] in the three-dimensional 
# space where the flow fields have been configured

(x,y,z)=(np*1/3,np*1/3,nz*1/3)

𝐼=Individuals(𝐹,x,y,z)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 1.4 A Closer Look (optional)
#
# The above `Individuals` constructor wraps up 𝐹, the initial position, and other needed components 
# within 𝐼. **At this point, you can either jump to section 2 or read through this section**
# to learn more about how the details as needed e.g. if you wanted to overide default options 
# that were selected for you by the section 1.3 constructor.
#
# Initial position is 

📌=[x,y,z] 

# and the data structure ([DataFrame](http://juliadata.github.io/DataFrames.jl/stable/)) 
# to record properties along the individual's path accordingly. 

🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], t=Float64[])

# It is the postprocessing function's responsibility to provide the record. It is thus 
# important that this intermediary (`postproc`) be consistent with the solver setup (`sol`) 
# and the expected record format (`🔴`).

function postproc(sol,𝐹::FlowFields;id=missing,𝑇=missing)
    df=postprocess_xy(sol,𝐹,id=id,𝑇=𝑇)
    #add third coordinate
    z=sol[3,:]
    df.z=z[:]
    return df
end

# The velocity function `🚄` relies only on flow fields obtained from
# `𝐹` (which is defined above) to interpolate velocity at the specified
# space-time position (e.g. those of individuals). 

🚄 = dxyz_dt

# Now that every thing needed to carry out the computation is in place, 
# we wrap up the problem configuration in a struct (`Individuals`) which 
# links to the initial positions, flow fields, etc. all that will be 
# necessary to compute trajectories over time (`∫!(𝐼,𝑇)`).

#assemble as a NamedTuple:
I=(position=📌,record=🔴,velocity=🚄,
postprocessing=postproc,parameters=𝐹)

#construct Individuals from NamedTuple:
𝐼=Individuals(I)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2 Trajectory Simulations
#
# The `∫!` function call below returns the final positions & updates `𝐼.📌` accordingly. It also records properties observed along the trajectory in `𝐼.🔴`. 
# Simple methods to visualize the individual trajectory (plot or movie) are provided at the end.

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 2.1 Compute Trajectories

𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 2.2 Visualize Trajectories
#
# - define `myplot` convenience function
# - single plot example using `myplot`
# - (generate animation using `myplot`)

#md p=dirname(pathof(IndividualDisplacements))
#md include(joinpath(p,"../examples/recipes_plots.jl"));
#md nt=length(𝐼.🔴.x)

#md myplot(i)=plot(𝐼.🔴.x[1:i],𝐼.🔴.y[1:i],𝐼.🔴.z[1:i],linewidth=2,arrow = 2,
#md     title="Solid body rotation / Spiral example",leg=false,
#md     xaxis="x",yaxis="y",zaxis="z",xlims=(0,np),ylims=(0,np));

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Single plot example:

#md plt=myplot(nt)
#md scatter!(plt,[📌[1]],[📌[2]],[📌[3]])
#md #scatter!(plt,[𝐼.🔴.x[end]],[𝐼.🔴.y[end]],[𝐼.🔴.z[end]])
#md scatter!(plt,[𝐼.📌[1]],[𝐼.📌[2]],[𝐼.📌[3]])

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Animation example:
#
# ```
# p=Int(ceil(nt/100))
# anim = @animate for i ∈ 1:p:nt
#     myplot(i)
# end
#
# pth=tempdir()*"/"
# gif(anim, pth*"SolidBodyRotation.gif", fps = 15)
# ```

