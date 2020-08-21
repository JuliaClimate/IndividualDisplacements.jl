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

using OrdinaryDiffEq, Plots, DataFrames
using IndividualDisplacements, MeshArrays

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.2  Gridded Domain

np,nz=16,4 #horizontal and vertical domain size
Γ=simple_periodic_domain(np)
γ=Γ["XC"].grid;

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.3 Velocity Fields

#Solid-body rotation around central location ...
i=Int(np/2+1)
u=-(Γ["YG"].-Γ["YG"][1][i,i])
v=(Γ["XG"].-Γ["XG"][1][i,i])

#... plus a convergent term to / from central location
d=-0.01
u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])

#Replicate u,v in vertical dimension
uu=MeshArray(γ,γ.ioPrec,nz)
[uu[k]=u[1] for k=1:nz]
vv=MeshArray(γ,γ.ioPrec,nz)
[vv[k]=v[1] for k=1:nz]

#Vertical velocity component w
w=fill(1.0,MeshArray(γ,γ.ioPrec,nz))

#store parameters in a NamedTuple
𝑃=(u0=uu, u1=uu, v0=vv, v1=vv,w0=0.0*w, w1=-0.01*w, 𝑇=[0,19.95*2*pi], ioSize=(np,np,nz));

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.4 Initial Position

xy=[np*1/3,np*1/3,nz*1/3]

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 Solve For Particle Trajectory

tr = DataFrame( ID=[], x=[], y=[], z=[], t = [])
solv(prob) = solve(prob,Tsit5(),reltol=1e-8)
function postproc(sol,𝑃::NamedTuple;id=missing,𝑇=missing)
    df=postprocess_xy(sol,𝑃,id=id,𝑇=𝑇)
    #add third coordinate
    z=sol[3,:,:]
    df.z=z[:]
    return df
 end

𝐼 = Individuals{Float64}(📌=xy[:,:], 🔴=tr, 🆔=collect(1:size(xy,2)),
                         🚄 = dxyz_dt, ∫ = solv, 🔧 = postproc, 𝑃=𝑃)
𝑇=(0.0,𝐼.𝑃.𝑇[2])
∫!(𝐼,𝑇)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.2 Visualize Particle Trajectory
#
# - define `myplot` convenience function
# - generate animation using `myplot`
# - single plot example using `myplot`

x,y,z=𝐼.🔴.x,𝐼.🔴.y,𝐼.🔴.z

myplot(i)=plot(x[1:i],y[1:i],z[1:i],linewidth=2,arrow = 2,
    title="Solid body rotation / Spiral example",leg=false,
    xaxis="x",yaxis="y",zaxis="z",xlims=(0,np),ylims=(0,np))

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Animation example:

nt=length(x)
p=Int(ceil(nt/100))
anim = @animate for i ∈ 1:p:nt
    myplot(i)
end
pth=tempdir()*"/"
gif(anim, pth*"SolidBodyRotation.gif", fps = 15)

#nb # %% {"slideshow": {"slide_type": "subslide"}}
# Single plot example:

plt=myplot(nt)
scatter!(plt,[xy[1]],[xy[2]],[xy[3]])
scatter!(plt,[x[end]],[y[end]],[z[end]])
