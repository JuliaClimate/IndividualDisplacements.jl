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

# Put grid variables in a dictionary.

# +
mygrid=gcmgrid("llc90_latlon/","ll",1,[(360,178)], [360 178], Float32, read, write)

GridVariables=Dict("XC" => read(mygrid.path*"XC.latlon.data",MeshArray(mygrid,Float32)),
"YC" => read(mygrid.path*"YC.latlon.data",MeshArray(mygrid,Float32)),
"XG" => read(mygrid.path*"XG.data",MeshArray(mygrid,Float32)),
"YG" => read(mygrid.path*"YG.data",MeshArray(mygrid,Float32)),
"DXC" => read(mygrid.path*"DXC.latlon.data",MeshArray(mygrid,Float32)),
"DYC" => read(mygrid.path*"DYC.latlon.data",MeshArray(mygrid,Float32)) );
# -

# Read velocity fields as `MeshArray`s.

# +
using MAT
file = matopen(mygrid.path*"uv_lonlat.mat")
u=read(file, "u")
v=read(file, "v")
close(file)

u=dropdims(mean(u,dims=3),dims=3)
v=dropdims(mean(v,dims=3),dims=3)

u=read(u,MeshArray(mygrid,Float32))
v=read(v,MeshArray(mygrid,Float32));

u[findall(isnan.(u))]=0.0
v[findall(isnan.(v))]=0.0

u0=u; u1=u;
v0=v; v1=v;
# -

# Put velocity fields and time range in a dictionary.

# +
t0=0.0; t1=86400*366*2.0; dt=3600;

u0=u0./GridVariables["DXC"]
u1=u1./GridVariables["DXC"]
v0=v0./GridVariables["DYC"]
v1=v1./GridVariables["DYC"]

uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1, "dt" => dt) ;
# -

# Merge the two dictionaries and add masks

# +
uvetc=merge(uvt,GridVariables);

nr=50; kk=1;

mskW=read(mygrid.path*"hFacW.latlon.data",MeshArray(mygrid,Float32,nr))
mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
mskS=read(mygrid.path*"hFacS.latlon.data",MeshArray(mygrid,Float32,nr))
mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)

msk=Dict("mskW" => mskW, "mskS" => mskS)

uvetc=merge(uvetc,msk);
# -

# ## 3. Compute trajectories from gridded flow fields

comp_vel=IndividualDisplacements.VelComp
get_vel=IndividualDisplacements.VelCopy

# +
ii1=0:5:360; ii2=20:2:150;

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

tspan = (0.0,uvt["t1"]-uvt["t0"])
prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
size(sol)

# ## 4. Store trajectories in a `DataFrame` and plot

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
lon=mod.(sol[1,:,:],360)
lat=mod.(sol[2,:,:],180)
df = DataFrame(ID=Int.(ID[:]), lon=lon[:], lat=lat[:])
size(df)

PyPlot.figure(); PlotBasic(df,size(sol,2),90.0)


