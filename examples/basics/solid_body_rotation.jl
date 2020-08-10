# # Single Particle Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/solid_body_rotation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/solid_body_rotation.ipynb)
#
# Simulate the trajectory of a particle in a perfectly circular flow (i.e.
# _solid body rotation_), which may represent e.g. an ocean meso-scale eddy.
# _Addendum: _ a homogeneous sinking / floating term was later added.
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
# - setup the software and initialize example
# - simulate trajectories & plot results
# - experiment with parameters (user)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software

using OrdinaryDiffEq, Plots
using IndividualDisplacements, MeshArrays

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.2  Gridded Domain
#
# - define `SetPeriodicDomain` function, which uses `MeshArrays.jl`
# - call `SetPeriodicDomain` function with a chosen grid size; e.g. `np=16` in
#   the horizontal directions and `nz=4` in the vertical.

np,nz=16,4
Γ=simple_periodic_domain(np);

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.3 Time & Velocity Fields
#

#time range
t0=0.0
t1=0.95*2*pi
#t1=2.95*2*pi
t1=19.95*2*pi

#solid-body rotation around central location
i=Int(np/2+1)
u=-(Γ["YG"].-Γ["YG"][1][i,i])
v=(Γ["XG"].-Γ["XG"][1][i,i])

#add some convergence to / divergence from central location
d=0.0
d=-0.01
u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])

#"vertical" component w
γ=Γ["XC"].grid
w=fill(1.0,MeshArray(γ,γ.ioPrec,nz))

#replicate u,v "vertically"
uu=MeshArray(γ,γ.ioPrec,nz)
[uu[k]=u[1] for k=1:nz]
vv=MeshArray(γ,γ.ioPrec,nz)
[vv[k]=v[1] for k=1:nz]

#store everything in a data structure
𝑃=(u0=uu, u1=uu, v0=vv, v1=vv,w0=0.0*w, w1=-0.01*w, 𝑇=[t0,t1], ioSize=(np,np,nz));

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.4 Initial Position and Time

u0=[np*1/3,np*1/3,nz*1/3]
du=fill(0.0,3)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Solve For Particle Trajectory
#
# - `ODEProblem` formulates the differential equation along with the time period `𝑇`, parameters `𝑃`
# - `solve` then performs the integration over `𝑇`, starting from `u0`
#
# _For additional documentation, try `?ODEProblem` or `?solve`_

prob = ODEProblem(dxyz_dt,u0,𝑃.𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-8)

x,y,z=sol[1,:],sol[2,:],sol[3,:]
nt=length(x)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.2 Visualize Particle Trajectory
#
# - define `myplot` convenience function
# - generate animation using `myplot`
# - single plot example using `myplot`

myplot(i)=plot(x[1:i],y[1:i],z[1:i],linewidth=2,arrow = 2,
    title="Solid body rotation / Spiral example",leg=false,
    xaxis="x",yaxis="y",zaxis="z",xlims=(0,np),ylims=(0,np))

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Animation example:

p=Int(ceil(nt/100))
anim = @animate for i ∈ 1:p:nt
    myplot(i)
end
pth=tempdir()*"/"
gif(anim, pth*"SolidBodyRotation.gif", fps = 15)

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Single plot example:

plt=myplot(nt)
scatter!(plt,[u0[1]],[u0[2]],[u0[3]])
scatter!(plt,[x[end]],[y[end]],[z[end]])
