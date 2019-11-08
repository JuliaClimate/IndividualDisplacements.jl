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
import IndividualDisplacements: myread

mygrid=gcmgrid("flt_example/","ll",1,[(80,42)], [80 42], Float32, read, write)

GridVariables=Dict("XC" => myread(mygrid.path*"XC",MeshArray(mygrid,Float32)),
"YC" => myread(mygrid.path*"YC",MeshArray(mygrid,Float32)),
"XG" => myread(mygrid.path*"XG",MeshArray(mygrid,Float32)),
"YG" => myread(mygrid.path*"YG",MeshArray(mygrid,Float32)),
"DXC" => fill(5000.0,MeshArray(mygrid,Float32)),
"DYC" => fill(5000.0,MeshArray(mygrid,Float32)));
# -

# Read velocity fields as `MeshArray`s.

# +
nr=8
u0=myread(mygrid.path*"U.0000000001",MeshArray(mygrid,Float32,nr))
u1=myread(mygrid.path*"U.0000018001",MeshArray(mygrid,Float32,nr))
v0=myread(mygrid.path*"V.0000000001",MeshArray(mygrid,Float32,nr))
v1=myread(mygrid.path*"V.0000018001",MeshArray(mygrid,Float32,nr))

kk=3 #3 to match -1406.25 in pkg/flt output
u0=u0[:,kk]; u1=u1[:,kk];
v0=v0[:,kk]; v1=v1[:,kk];
# -

# Put velocity fields and time range in a dictionary.

# +
t0=0.0 #approximation / simplification
t1=18001.0*3600.0
dt=1

u0=u0./GridVariables["DXC"]
u1=u1./GridVariables["DXC"]
v0=v0./GridVariables["DYC"]
v1=v1./GridVariables["DYC"]

uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;
# -

# Merge the two dictionaries and add masks

# +
uvetc=merge(uvt,GridVariables);

mskW=myread(mygrid.path*"hFacW",MeshArray(mygrid,Float32,nr))
mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
mskS=myread(mygrid.path*"hFacS",MeshArray(mygrid,Float32,nr))
mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)

msk=Dict("mskW" => mskW, "mskS" => mskS)

uvetc=merge(uvetc,msk);
# -

# ## 3. Visualize  gridded variables

heatmap(mskW[1,1].*u0[1,1],title="U at the start")

# ## 4. Recompute displacements from gridded flow fields

# Initialize individual locations and define method aliases.

# +
uInit=[37.525990625,4.22379453125]
nSteps=2998
du=fill(0.0,2);

comp_vel=IndividualDisplacements.VelComp
get_vel=IndividualDisplacements.VelCopy
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
#tmp=zeros(10,1)
#for comp_vel(du,[180.1;40.1],uvetc,0.0)
Plots.plot(tmpu)
Plots.plot!(tmpv)

# ## 5. Solve through time using `DifferentialEquations.jl`

using DifferentialEquations
tspan = (0.0,nSteps*3600.0)
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
#ii1=1:10:80; ii2=1:10:42; #->sol is (2, 40, 40065)
#ii1=30:37; ii2=16:20; #->sol is (2, 40, 9674)
#ii1=10:17; ii2=16:20; #->sol is (2, 40, 51709)
ii1=5:5:40; ii2=5:5:25; #->sol is (2, 40, 51709)

n1=length(ii1); n2=length(ii2);
uInitS=Array{Float64,2}(undef,(2,n1*n2))
for i1 in eachindex(ii1); for i2 in eachindex(ii2);
        i=i1+(i2-1)*n1
        uInitS[1,i]=ii1[i1]-0.5
        uInitS[2,i]=ii2[i2]-0.5       
end; end;
du=fill(0.0,size(uInitS));
comp_vel(du,uInitS,uvetc,0.0)
du
# -

prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
size(sol)

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=5000* mod.(sol[1,:,:],80)
lat=5000* mod.(sol[2,:,:],42)
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:])
size(df)

PyPlot.figure(); PlotBasic(df,size(sol,2),100000.0)

Plots.plot(sol_one[1,:],sol_one[2,:])
#i1=38.0; i2=5.0; i=i1+(i2-1)*n1;
i=5
Plots.plot!(mod.(sol[1,i,:],80),mod.(sol[2,i,:],42))


