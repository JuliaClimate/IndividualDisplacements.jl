# # Particle Cloud
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/particle_cloud.ipynb)
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

using IndividualDisplacements, Statistics

import IndividualDisplacements.OrdinaryDiffEq: solve, Tsit5
import IndividualDisplacements.DataFrames: DataFrame

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/jupyter/example123.jl"));
#md include(joinpath(p,"../examples/jupyter/recipes_plots.jl"))

# ## 2. Setup Problem

𝑃,Γ=example2_setup()

ii1=5:5:40; ii2=5:5:25
x=vec([x-0.5 for x in ii1, y in ii2])
y=vec([y-0.5 for x in ii1, y in ii2])
xy = permutedims([[x[i];y[i];1.0] for i in eachindex(x)])

solv(prob) = IndividualDisplacements.ensemble_solver(prob,solver=Tsit5(),reltol=1e-6,abstol=1e-6)
tr = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])

#𝐼 = Individuals{Float64,2}(📌=xy[:,:], 🔴=tr, 🆔=collect(1:size(xy,2)),
#                         🚄 = dxdt!, ∫ = solv, 🔧 = postprocess_xy, 𝑃=𝑃);

I=(position=xy,record=deepcopy(tr),velocity=dxdt!,
   integration=solv,postprocessing=postprocess_xy,parameters=𝑃)
𝐼=Individuals(I)

# ## 3. Compute Trajectories

𝑇 = (0.0,2998*3600.0)
∫!(𝐼,𝑇)

# ## 4. Display results

#md 𝐼.🔴.lon=5000*𝐼.🔴.x
#md 𝐼.🔴.lat=5000*𝐼.🔴.y
#md plt=plot_paths(𝐼.🔴,size(xy,2),100000.0)

# Compare with trajectory output from `MITgcm`

#IndividualDisplacements.flt_example_download()
#df=read_flt(IndividualDisplacements.flt_example_path,Float32)
#ref=plot_paths(df,size(xy,2),100000.0)
