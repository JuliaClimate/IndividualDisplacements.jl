# # Detailed Look
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/detailed_look.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/detailed_look.ipynb)
#
# A more detailed look at spatial interpolation, integration through time, and I/O. 
# For additional documentation e.g. see
# [1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector)). 
# Here we illustrate a few things in more detail:
#
# - reading velocities from file.
#   - gridded velocity output (U*data, V*data)
#   - pre-computed trajectory output (`float_traj*data`)
# - interpolating `U,V` from gridded output to individual locations
#   - compared with `u,v` from `float_traj*data`
# - computing trajectories (location v time) using `OrdinaryDiffEq.jl`
#   - compared with `x(t),y(t)` from `float_traj*data`

# ## 1. Import Software

using IndividualDisplacements, OrdinaryDiffEq, DataFrames, MITgcmTools
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))
include(joinpath(p,"../examples/example123.jl"))
include(joinpath(p,"../examples/helper_functions.jl"));

# ## 2. Read Trajectory Output
#
# from `MITgcm/pkg/flt`

dirIn=IndividualDisplacements.flt_example
prec=Float32
df=read_flt(dirIn,prec);

#md plt=plot_paths(df,300,100000.0)

# ## 3. Read Gridded Variables
#
# using `MeshArrays.jl` and e.g. a NamedTyple

𝑃,Γ=example2_setup();

# ## 4. Visualize Velocity Fields
#
# ```
# plt=heatmap(Γ["mskW"][1,1].*𝑃.u0,title="U at the start")
# plt=heatmap(Γ["mskW"][1,1].*𝑃.u1-𝑃.u0,title="U end - U start")
# ```

# ## 5. Visualize Trajectories
#
# (select one trajectory)

tmp=df[df.ID .== 200, :]
tmp[1:4,:]

# Super-impose trajectory over velocity field (first for u ...)

x=Γ["XG"].f[1][:,1]
y=Γ["YC"].f[1][1,:]
z=transpose(Γ["mskW"][1].*𝑃.u0);

# plt=contourf(x,y,z,c=:delta)
# plot!(tmp[:,:lon],tmp[:,:lat],c=:red,w=4,leg=false)

# Super-impose trajectory over velocity field (... then for v)

x=Γ["XC"].f[1][:,1]
y=Γ["YG"].f[1][1,:]
z=transpose(Γ["mskW"][1].*𝑃.v0);

# plt=contourf(x,y,z,c=:delta)
# plot!(tmp[:,:lon],tmp[:,:lat],c=:red,w=4,leg=false)

# ## 6. Interpolate Velocities

dx=Γ["dx"]
uInit=[tmp[1,:lon];tmp[1,:lat]]./dx
nSteps=Int32(tmp[end,:time]/3600)-2
du=fill(0.0,2);

# Visualize and compare with actual grid point values -- jumps on the tangential component are expected with linear scheme:

tmpu=fill(0.0,100)
tmpv=fill(0.0,100)
tmpx=fill(0.0,100)
for i=1:100
    tmpx[i]=500.0 *i./dx
    dxy_dt(du,[tmpx[i];0.499./dx],𝑃,0.0)
    tmpu[i]=du[1]
    tmpv[i]=du[2]
end

#md plt=plot(tmpx,tmpu,label="u (interp)")
#md plot!(Γ["XG"].f[1][1:10,1]./dx,𝑃.u0[1:10,1],marker=:o,label="u (C-grid)")
#md plot!(tmpx,tmpv,label="v (interp)")
#md plot!(Γ["XG"].f[1][1:10,1]./dx,𝑃.v0[1:10,1],marker=:o,label="v (C-grid)")

# And similarly in the other direction

tmpu=fill(0.0,100)
tmpv=fill(0.0,100)
tmpy=fill(0.0,100)
for i=1:100
    tmpy[i]=500.0 *i./dx
    dxy_dt(du,[0.499./dx;tmpy[i]],𝑃,0.0)
    tmpu[i]=du[1]
    tmpv[i]=du[2]
end

# plt=plot(tmpx,tmpu,label="u (interp)")
# plot!(Γ["YG"].f[1][1,1:10]./dx,𝑃.u0[1,1:10],marker=:o,label="u (C-grid)")
# plot!(tmpx,tmpv,label="v (interp)")
# plot!(Γ["YG"].f[1][1,1:10]./dx,𝑃.v0[1,1:10],marker=:o,label="v (C-grid)")

# Compare recomputed velocities with those from `pkg/flt`

nSteps=2998
tmpu=fill(0.0,nSteps); tmpv=fill(0.0,nSteps);
tmpx=fill(0.0,nSteps); tmpy=fill(0.0,nSteps);
refu=fill(0.0,nSteps); refv=fill(0.0,nSteps);
for i=1:nSteps
    dxy_dt_replay(du,[tmp[i,:lon],tmp[i,:lat]],tmp,tmp[i,:time])
    refu[i]=du[1]./dx
    refv[i]=du[2]./dx
    dxy_dt(du,[tmp[i,:lon],tmp[i,:lat]]./dx,𝑃,Float64(tmp[i,:time]))
    tmpu[i]=du[1]
    tmpv[i]=du[2]
end

# plt=plot(tmpu,label="u")
# plot!(tmpv,label="v")
# plot!(refu,label="u (ref)")
# plot!(refv,label="v (ref)")

# ## 6. Compute Trajectories
#
# Solve through time using `OrdinaryDiffEq.jl` with
#
# - `dxy_dt` is the function computing `dxy/dt`
# - `uInit` is the initial condition `u @ tspan[1]`
# - `tspan` is the time interval
# - `uvetc` are parameters for `dxy_dt`
# - `Tsit5` is the time-stepping scheme
# - `reltol` and `abstol` are tolerance parameters

tspan = (0.0,nSteps*3600.0)
#prob = ODEProblem(dxy_dt_replay,uInit,tspan,tmp)
prob = ODEProblem(dxy_dt,uInit,tspan,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
sol[1:4]

# Compare recomputed trajectories with originals from `MITgcm/pkg/flt`

ref=transpose([tmp[1:nSteps,:lon] tmp[1:nSteps,:lat]])
maxLon=80*5.e3
maxLat=42*5.e3
show(size(ref))
for i=1:nSteps-1
    ref[1,i+1]-ref[1,i]>maxLon/2 ? ref[1,i+1:end]-=fill(maxLon,(nSteps-i)) : nothing
    ref[1,i+1]-ref[1,i]<-maxLon/2 ? ref[1,i+1:end]+=fill(maxLon,(nSteps-i)) : nothing
    ref[2,i+1]-ref[2,i]>maxLat/2 ? ref[2,i+1:end]-=fill(maxLat,(nSteps-i)) : nothing
    ref[2,i+1]-ref[2,i]<-maxLat/2 ? ref[2,i+1:end]+=fill(maxLat,(nSteps-i)) : nothing
end
ref=ref./dx;

#md plt=plot(sol[1,:],sol[2,:],linewidth=5,title="Using Recomputed Velocities",
#md      xaxis="lon",yaxis="lat",label="Julia Solution") # legend=false
#md plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,label="MITgcm Solution")
