# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# # This notebook
#
# _Notes:_ For documentation see <https://gaelforget.github.io/MeshArrays.jl/stable/>, <https://docs.juliadiffeq.org/latest/solvers/ode_solve.html> and <https://en.wikipedia.org/wiki/Displacement_(vector)>

# ## 1. import software

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Plots, Statistics, MITgcmTools, DataFrames
include(joinpath(dirname(pathof(MeshArrays)),"../examples/Plots.jl"))

# ## 2. Read gridded variables as `MeshArray`s

# Put grid variables in a dictionary.

γ=GridSpec("LatLonCap","GRID_LLC90/")
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))
uvetc=read_uv_etc(20,Γ)

# Visualize  gridded variables

plt=heatmap(uvetc["u0"][1,1],title="U at the start")
display(plt)

# ## 3. Compute trajectories from gridded flow fields

# Set `comp_vel`, an alias, to a suitable function.

comp_vel=IndividualDisplacements.VelComp!

# Inspect how `comp_vel` behaves in small test case.

# +
uInit=[45.0,100.0,1.0]
du=fill(0.0,3);

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
plt=Plots.plot(tmpu)
Plots.plot!(tmpv)
display(plt)
# -

# Solve for trajectory in small test case.

tspan = (0.0,uvt["t1"]-uvt["t0"])
prob = ODEProblem(comp_vel,uInit,tspan,uvetc)
sol_one = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4)
sol_two = solve(prob,Euler(),dt=1e6)
size(sol_one)

# Define initial condition array

# +
uInitS = Array{Float64,2}(undef, 3, prod(XC.grid.ioSize))
kk = 0

for fIndex = 1:5
        nx, ny = XC.fSize[fIndex]
        ii1 = 0.5:1.0:nx
        ii2 = 0.5:1.0:ny
        n1 = length(ii1)
        n2 = length(ii2)
        for i1 in eachindex(ii1)
          for i2 in eachindex(ii2)
            if msk[fIndex][Int(round(i1+0.5)),Int(round(i2+0.5))]
                        global kk += 1
                        let kk = kk
                                uInitS[1, kk] = ii1[i1]
                                uInitS[2, kk] = ii2[i2]
                                uInitS[3, kk] = fIndex
                        end
            end
          end
        end
end

uInitS=uInitS[:,1:kk]
du=fill(0.0,size(uInitS));
# -

# Solve for all trajectories.

prob = ODEProblem(comp_vel,uInitS,tspan,uvetc)
sol = solve(prob,Euler(),dt=uvt["dt"])
size(sol)

# ## 4. Plot trajectories
#
# - Copy `sol` to a `DataFrame`

ID=collect(1:size(sol,2))*ones(1,size(sol,3))
x=sol[1,:,:]
y=sol[2,:,:]
fIndex=sol[3,:,:]
df = DataFrame(ID=Int.(ID[:]), x=x[:], y=y[:], fIndex=fIndex[:])
size(df)

# - Map i,j position to lon,lat coordinates

# +
lon=Array{Float64,1}(undef,size(df,1))
lat=similar(lon)

for ii=1:length(lon)

#get location in grid index space
x=df[ii,:x]; y=df[ii,:y]; fIndex=Int(df[ii,:fIndex])
dx,dy=[x - floor(x),y - floor(y)]
i_c,j_c = Int32.(floor.([x y])) .+ 2

#interpolate lon and lat to position
tmp=view(YC[fIndex],i_c:i_c+1,j_c:j_c+1)
lat[ii]=(1.0-dx)*(1.0-dy)*tmp[1,1]+dx*(1.0-dy)*tmp[2,1]+
(1.0-dx)*dy*tmp[1,2]+dx*dy*tmp[2,2]

tmp=view(XC[fIndex],i_c:i_c+1,j_c:j_c+1)
if (maximum(tmp)>minimum(tmp)+180)&&(lat[ii]<88)
        tmp1=deepcopy(tmp)
        tmp1[findall(tmp.<maximum(tmp)-180)] .+= 360.
        tmp=tmp1
end
#kk=findall(tmp.<maximum(tmp)-180); tmp[kk].=tmp[kk].+360.0
lon[ii]=(1.0-dx)*(1.0-dy)*tmp[1,1]+dx*(1.0-dy)*tmp[2,1]+
(1.0-dx)*dy*tmp[1,2]+dx*dy*tmp[2,2]
end

df.lon=lon
df.lat=lat

show(df)
# -

# - call `PlotMapProj`

# +
p=dirname(pathof(IndividualDisplacements))

nn=10000

#include(joinpath(p,"plot_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

include(joinpath(p,"plot_makie.jl"))
AbstractPlotting.inline!(true) #for Juno, set to false
scene=PlotMakie(df,nn,180.0)
#Makie.save("LatLonCap300mDepth.png", scene)
