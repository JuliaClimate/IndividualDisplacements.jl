# # Single Particle Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/solid_body_rotation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/solid_body_rotation.ipynb)
#
# Simulate the trajectory of a particle in a perfectly circular flow (i.e.
# _solid body rotation_), which may represent e.g. an ocean meso-scale eddy.
#
# ![solid body rotation](https://github.com/JuliaClimate/IndividualDisplacements.jl/raw/master/examples/figs/SolidBodyRotation.gif)
#
# As an exercise left to the user, directions are provided e.g. to add a convergence / divergence term.
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# - 1. setup the software and initialize example
# - 2. simulate trajectories & plot results
# - 3. experiment with parameters (user)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software

using OrdinaryDiffEq, Plots
using IndividualDisplacements, MeshArrays

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.2  Gridded Domain
#
# - define `SetPeriodicDomain` function, which uses `MeshArrays.jl`
# - call `SetPeriodicDomain` function with a chosen grid size; e.g. `np=16`

np=16

Γ=simple_periodic_domain(np);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.3 Time & Velocity Fields
#
# - define time range
# - define velocity field(s)
# - store in `𝑃` (dictionary) with grid variables

#time range
t0=0.0
t1=0.95*2*pi
#t1=2.95*2*pi

#solid-body rotation around central location
i=Int(np/2+1)
u=-(Γ["YG"].-Γ["YG"][1][i,i])
v=(Γ["XG"].-Γ["XG"][1][i,i])

#add some convergence to / divergence from central location
d=0.0
#d=-0.10
u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])

#store everything in a dictionnary
𝑃=Dict("u0" => u, "u1" => u, "v0" => v, "v1" => v, "t0" => t0, "t1" => t1)
𝑃=merge(𝑃,Γ);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.4 Initial Position and Time

u0=np*[1/3,1/3]
du=fill(0.0,2)
𝑇 = (𝑃["t0"],𝑃["t1"]);

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Solve For Particle Trajectory
#
# - `ODEProblem` formulates the differential equation along with the time period `𝑇`, parameters `𝑃`
# - `solve` then performs the integration over `𝑇`, starting from `u0`
#
# _For additional documentation, try `?ODEProblem` or `?solve`_

prob = ODEProblem(⬡,u0,𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-8)

x,y=sol[1,:],sol[2,:]
nt=length(x)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.2 Visualize Particle Trajectory
#
# - define `myplot` convenience function
# - generate animation using `myplot`
# - single plot example using `myplot`

myplot(i)=plot(x[1:i],y[1:i],linewidth=2,arrow = 2,
    title="Solid body rotation / Spiral example",leg=false,
    xaxis="x",yaxis="y",xlims=(0,np),ylims=(0,np))

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Animation example:

if false
p=Int(ceil(nt/100))
anim = @animate for i ∈ 1:p:nt
    myplot(i)
end
pth=tempdir()*"/"
gif(anim, pth*"SolidBodyRotation.gif", fps = 15)
end

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Single plot example:

plt=myplot(nt)
scatter!(plt,[u0[1]],[u0[2]])
