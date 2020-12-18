# # Single Particle
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
# ![solid body rotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# # 1 Problem Configuration
#
# Here we set up software, grid, flow fields, initial conditions.
#
# ### 1.1 Import Software

using IndividualDisplacements, DataFrames
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### 1.2  Gridded Domain

np,nz=16,4 #horizontal and vertical domain size
Γ=simple_periodic_domain(np);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### 1.3 Velocity Fields
#
# Exercise: find `simple_flow_field` within `helper_functions.jl` and modify the 
# flow field parameters (e.g. intensity and sign of the convergent term).

u,v,w=simple_flow_field(Γ,np,nz);

# ### 1.4 Velocity Function
#
# `🚄` relies only on parameters (velocity fields, grid, etc) 
# contained in `𝑄` to compute velocity at the space-time position
# of the individual. The solver (here: `solv`) can then integrate 
# over time the result of `🚄` (see `OrdinaryDiffEq.jl` docs).

🚄 = dxyz_dt

𝑄=𝑃_Array3D{eltype(u)}(u,u,v,v,0*w,1*w,[0,19.95*2*pi])

solv(prob) = solve(prob,Tsit5(),reltol=1e-8)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 1.5 Initial Position
#
# Here we set up just one individual in a three-dimensional space,

📌=[np*1/3,np*1/3,nz*1/3]

# and the data structure ([DataFrame](http://juliadata.github.io/DataFrames.jl/stable/)) 
# to record properties along the individual's path accordingly. It is the postprocessing 
# function's responsibility to provide the record. It is thus important that this 
# intermediary (`postproc`) be consistent with the solver setup (`sol`) and 
# the expected record format (`🔴`).

🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], t=Float64[])

function postproc(sol,𝑄::FlowParameters;id=missing,𝑇=missing)
    df=postprocess_xy(sol,𝑄,id=id,𝑇=𝑇)
    #add third coordinate
    z=sol[3,:]
    df.z=z[:]
    return df
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2 Trajectory Simulations
#
# Now that every thing needed to carry out the computation is in place, 
# we wrap up the problem configuration in a struct (`Individuals`) which 
# links to the initial positions, flow fields, etc. all that will be 
# necessary to compute trajectories over time (`∫!(𝐼,𝑇)`). Simple methods to
# visualize the individual trajectory (plot or movie) are provided at the end.

# ### 2.1 Setup Individuals
#

#assemble as a NamedTuple:
I=(position=📌,record=🔴,velocity=🚄,
integration=solv,postprocessing=postproc,parameters=𝑄)

#construct Individuals from NamedTuple:
𝐼=Individuals(I)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 2.2 Compute Trajectories
#
# The `∫!` function call below returns the final positions & updates `𝐼.📌` accordingly. It also records properties observed along the trajectory in `𝐼.🔴`

𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### 2.3 Visualize Trajectories
#
# - define `myplot` convenience function
# - generate animation using `myplot`
# - single plot example using `myplot`

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

#md p=Int(ceil(nt/100))
#md anim = @animate for i ∈ 1:p:nt
#md     myplot(i)
#md end

#md pth=tempdir()*"/"
#md gif(anim, pth*"SolidBodyRotation.gif", fps = 15)

# Exercise: make the sinking velocity decrease with time 
# (hint: it increases as specified above in the original notebook); 
# change the number of times the particle goes around the origin; etc
