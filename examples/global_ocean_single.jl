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
t0=0.0; t1=86400*366*10.0; dt=3600;

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

# ## 3. Visualize  gridded variables

#PyPlot.contourf(GridVariables["XG"].f[1], GridVariables["YC"].f[1], mskS.f[1].*v0.f[1])
#colorbar()
heatmap(v0.f[1])
#heatmap(mskS.f[1].*v0.f[1])
#PyPlot.contour(GridVariables["XG"].f[1], GridVariables["YC"].f[1], mskS.f[1].*v0.f[1])

# ## 4. Recompute displacements from gridded flow fields

# Initialize individual locations and define method aliases.

# +
uInit=[180.0,40.0] #uInit=[tmp[1,:lon];tmp[1,:lat]]./uvetc["dx"]
nSteps=Int32(uvt["t1"]/uvt["dt"]) #nSteps=Int32(tmp[end,:time]/3600)-2
du=fill(0.0,2);

comp_vel=IndividualDisplacements.VelComp
get_vel=IndividualDisplacements.VelCopy
# -

# Inspect how `comp_vel` behaves.

# +
if false 
    comp_vel(du,uInit,uvetc,0.0); show(du)
    tmpdu=[uvt["u0"][1][uInit[1]+1,uInit[2]+1] uvt["v0"][1][uInit[1]+1,uInit[2]+1]]; show(tmpdu)
end

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
# -

# ## 5. Solve through time using `DifferentialEquations.jl`

using DifferentialEquations
tspan = (0.0,nSteps*3600.0)
#prob = ODEProblem(get_vel,uInit,tspan,tmp)
prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
sol[1:4]

?Tsit5

#Plots.plot(sol[1,:],sol[2,:])
sol[:,end-4:end]

Plots.plot(sol[1,2:end]-sol[1,1:end-1])
Plots.plot!(sol[2,2:end]-sol[2,1:end-1])

Plots.plot(sol[1,:],sol[2,:])

