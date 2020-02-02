# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.0-rc4
#     language: julia
#     name: julia-1.3
# ---

# # This notebook
#
# _Notes:_ For documentation see <https://gaelforget.github.io/MeshArrays.jl/stable/>, <https://docs.juliadiffeq.org/latest/solvers/ode_solve.html> and <https://en.wikipedia.org/wiki/Displacement_(vector)>

# ## 1. import software

using IndividualDisplacements, MeshArrays, DifferentialEquations
using DataFrames, Plots, Statistics
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))

# ## 2. Define gridded variables as `MeshArray`s

# Put grid variables in a dictionary.

GridVariables=GridOfOnes("PeriodicDomain",1,40)
(Rini,Rend,DXCsm,DYCsm)=demo2(GridVariables)

heatmap(Rend[1])

# Derive velocity fields using Rend as a scalar potential (**later: streamfunction...**)

(u,v)=gradient(Rend,GridVariables)
u0=-v; u1=-v; v0=u; v1=u;

# Put velocity fields and time range in a dictionary.

# +
t0=0.0 #approximation / simplification
t1=100.0
dt=0.1
nSteps=(t1-t0)/dt

u0=u0./GridVariables["DXC"]#normalization to grid units
u1=u1./GridVariables["DXC"]
v0=v0./GridVariables["DYC"]
v1=v1./GridVariables["DYC"]

uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;
# -

# Merge the two dictionaries and add masks

# +
uvetc=merge(uvt,GridVariables);

mskW=fill(1.0,u)
mskS=fill(1.0,v)

msk=Dict("mskW" => mskW, "mskS" => mskS)

uvetc=merge(uvetc,msk);
# -

# ## 3. Visualize  gridded variables

heatmap(mskW[1,1].*u0[1,1],title="U at the start")

# ## 4. Recompute displacements from gridded flow fields

# Initialize individual locations and define method aliases.

# +
uInit=[20.0,20.0]
du=fill(0.0,2);

comp_vel=IndividualDisplacements.VelComp
# -

# Inspect how `comp_vel` behaves.

ii=uInit[1]-3:0.1:uInit[1]+3
jj=uInit[2]-3:0.1:uInit[2]+3
tmpu=zeros(size(ii))
tmpv=zeros(size(ii))
for i in eachindex(ii)
    comp_vel(du,[ii[i];jj[i]],uvetc,0.0)
    tmpu[i],tmpv[i]=du
end
Plots.plot(tmpu)
Plots.plot!(tmpv)

# ## 5. Solve through time using `DifferentialEquations.jl`

tspan = (0.0,nSteps*dt)
prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
size(sol_one)

uInitS=[uInit uInit]
du=fill(0.0,size(uInitS));
comp_vel(du,uInitS,uvetc,0.0)
du

# +
ii1=0.25:0.25:40; ii2=0.25:0.25:40;

n1=length(ii1); n2=length(ii2);
uInitS=Array{Float64,2}(undef,(2,n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2);
        i=i1+(i2-1)*n1
        uInitS[1,i]=ii1[i1]
        uInitS[2,i]=ii2[i2]
end; end;
du=fill(0.0,size(uInitS));
comp_vel(du,uInitS,uvetc,0.0)
du
# -

prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
size(sol)

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=mod.(sol[1,:,:],40)
lat=mod.(sol[2,:,:],40)
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:])
size(df)

# +
nn=minimum([5000 size(du,2)])

#p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"plot_pyplot.jl"))
#PyPlot.figure(); PlotBasic(df,nn,20.0)

p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"plot_makie.jl"))
AbstractPlotting.inline!(true) #for Juno, set to false
scene=PlotMakie(df,nn,20.0)
#Makie.save("PeriodicDomainRandomFlow_Makie.png", scene)
