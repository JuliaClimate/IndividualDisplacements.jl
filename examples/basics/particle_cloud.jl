# # Particle Cloud Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/particle_cloud.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/particle_cloud.ipynb)
#
# Using the same setup as `detailed_look.jl` or `example2()`, here we simulate
# a point cloud getting advected by the flow field. 
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))

# ## 1. Import Software

using IndividualDisplacements, OrdinaryDiffEq, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))
include(joinpath(p,"../examples/example123.jl"));

# ## 2. Setup Problem

𝑃,Γ=example2_setup()

ii1=5:5:40; ii2=5:5:25
x=vec([x-0.5 for x in ii1, y in ii2])
y=vec([y-0.5 for x in ii1, y in ii2])
xy=transpose([x y])

𝑃.𝑇[:] = [0.0,2998*3600.0]
solv(prob) = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6)
tr = DataFrame( ID=[], x=[], y=[], t = [])

𝐼 = Individuals{Float64}(xy=xy[:,:], 𝑃=𝑃, ⎔ = dxy_dt, □ = solv, ▽ = postprocess_xy, tr = tr);

# ## 3. Compute Trajectories

start!(𝐼)

# ## 4. Display results

𝐼.tr.lon=5000*𝐼.tr.x
𝐼.tr.lat=5000*𝐼.tr.y
plt=PlotBasic(𝐼.tr,size(xy,2),100000.0)

# Compare with trajectory output from `MITgcm`

#df=read_flt(joinpath(p,"../examples/flt_example/"),Float32)
#ref=PlotBasic(df,size(xy,2),100000.0)
