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

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq, Statistics, DataFrames
p=dirname(pathof(IndividualDisplacements)); include(joinpath(p,"plot_plots.jl"))
p=dirname(pathof(MeshArrays)); include(joinpath(p,"../examples/Demos.jl"))

# ## 2. Define gridded variables as `MeshArray`s

# Put grid variables in a dictionary.

# +
γ,Γ=GridOfOnes("PeriodicDomain",16,20)
(Rini,Rend,DXCsm,DYCsm)=demo2(Γ)

plt=heatmap(Rend[1])
display(plt)
# -

# Derive velocity fields using Rend as a scalar potential (**later: streamfunction...**)

(u,v)=gradient(Rend,Γ)
u=u./Γ["DXC"]#normalization to grid units
v=v./Γ["DYC"]
(u,v)=exchange(u,v,1)
u0=-v; u1=-v; v0=u; v1=u;

# Put velocity fields and time range in a dictionary.

# +
t0=0.0 #approximation / simplification
t1=200.0
dt=0.1
nSteps=(t1-t0)/dt

uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;
# -

# Merge the two dictionaries and add masks

# +
uvetc=merge(uvt,Γ);

# _Note:_ in general case, mskW & mskS would need to be exchanged ...

mskW=fill(1.0,u)
mskS=fill(1.0,v)

msk=Dict("mskW" => mskW, "mskS" => mskS)

uvetc=merge(uvetc,msk);
# -

# ## 3. Visualize  gridded variables

plt=heatmap(mskW[1,1].*u0[1,1],title="U at the start")
display(plt)

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
plot(tmpu)
plot!(tmpv)
# -

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
ii1=0.5:0.5:20; ii2=0.5:0.5:20; ii3=[2.0 4.0 5.0 7.0 10.0 12.0 13.0 15.0];

n1=length(ii1); n2=length(ii2); n3=length(ii3);
uInitS=Array{Float64,2}(undef,(3,n3*n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2); for i3 in eachindex(ii3);
        i=i1+(i2-1)*n1+(i3-1)*n1*n2
        uInitS[1,i]=ii1[i1]
        uInitS[2,i]=ii2[i2]
        uInitS[3,i]=ii3[i3]
end; end; end;
du=fill(0.0,size(uInitS));
comp_vel(du,uInitS,uvetc,0.0)
du
# -

prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
size(sol)

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=sol[1,:,:]+20.0*mod.(sol[3,:,:].-1,4)
lat=sol[2,:,:]+20.0*floor.((sol[3,:,:].-1)/4)
fIndex=sol[3,:,:]
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:], fIndex=fIndex[:])
size(df)

nn=minimum([5000 size(du,2)])
plt=PlotBasic(df,nn,10.0)
#savefig("PeriodicDomainRandomFlow.png")
display(plt)
