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
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# # This notebook
#
# _Notes:_ For documentation see <https://gaelforget.github.io/MeshArrays.jl/stable/>, <https://docs.juliadiffeq.org/latest/solvers/ode_solve.html> and <https://en.wikipedia.org/wiki/Displacement_(vector)>

# ## 1. import software

using IndividualDisplacements, MeshArrays, DifferentialEquations, Plots, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"PlotIndDisp.jl"))

# ## 2. Read gridded variables as `MeshArray`s

# _Note:_ `myread` function deals with tiled files from `flt_example/`; should be able to use it via `gcmgrid`...
#
# Put grid variables in a dictionary.

# +
#GridVariables=GridOfOnes("dpdo",1,40)
#(Rini,Rend,DXCsm,DYCsm)=MeshArrays.demo2(GridVariables)

GridVariables=GridOfOnes("ll",1,40)
(Rini,Rend,DXCsm,DYCsm)=MeshArrays.demo2(GridVariables)
heatmap(Rend[1])
# -

# Derive velocity fields using Rend as a scalar potential (**later: streamfunction...**)

(u,v)=gradient(Rend,GridVariables)
(u,v)=exchange(u,v,1)
u0=-v; u1=-v; v0=u; v1=u;

# +
dxc,dyc=exchange(GridVariables["DXC"],GridVariables["DYC"])
dxc=abs.(dxc)
dyc=abs.(dyc)

#hack to fix ll grid case
dxc[findall(dxc.<1.0)]=1.0
dyc[findall(dyc.<1.0)]=1.0

u0=u0./dxc#normalization to grid units
u1=u1./dxc
v0=v0./dyc
v1=v1./dyc
# -

# Put velocity fields and time range in a dictionary.

# +
t0=0.0 #approximation / simplification
t1=100.0
dt=0.1
nSteps=(t1-t0)/dt

uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;
# -

# Merge the two dictionaries and add masks

# +
uvetc=merge(uvt,GridVariables);

#this would need to be exchanged too...
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
uInit=[20.0,20.0,1.0]
du=fill(0.0,3);

comp_vel=IndividualDisplacements.VelComp!
get_vel=IndividualDisplacements.VelCopy
# -

# Inspect how `comp_vel` behaves.

# +
ii=uInit[1]-3:0.1:uInit[1]+3
jj=uInit[2]-3:0.1:uInit[2]+3
fIndex=ones(size(jj))

tmpu=zeros(size(ii))
tmpv=zeros(size(ii))
tmpf=zeros(size(ii))
for i in eachindex(ii)
    comp_vel(du,[ii[i];jj[i];fIndex[i]],uvetc,0.0)
    tmpu[i],tmpv[i],tmpf[i]=du
end
#tmp=zeros(10,1)
#for comp_vel(du,[180.1;40.1],uvetc,0.0)
Plots.plot(tmpu)
Plots.plot!(tmpv)
# -

# ## 5. Solve through time using `DifferentialEquations.jl`

using DifferentialEquations
tspan = (0.0,nSteps*dt)
#prob = ODEProblem(get_vel,uInit,tspan,tmp)
prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
size(sol_one)

uInitS=[uInit uInit]
du=fill(0.0,size(uInitS));
comp_vel(du,uInitS,uvetc,0.0)
du

# **Note how the size of sol (ie nb of steps) depends on initial location:**

# +
#ii1=1:40; ii2=1:40;
ii1=0.25:0.25:40; ii2=0.25:0.25:40;

n1=length(ii1); n2=length(ii2);
uInitS=Array{Float64,2}(undef,(3,n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2);
        i=i1+(i2-1)*n1
        uInitS[1,i]=ii1[i1]-0.5
        uInitS[2,i]=ii2[i2]-0.5       
        uInitS[3,i]=1.0
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

nn=minimum([5000 size(du,2)])
PyPlot.figure(); PlotBasic(df,nn,20.0)


