
module PlottingFunctions 

using GLMakie, DataFrames, FileIO, Colors
using DataDeps, MeshArrays

lon180(x)=Float64(x>180.0 ? x-360.0 : x)
lon360(x)=Float64(x<0.0 ? x+360.0 : x)

function background()
    dx=0.1
    lon,lat,basemap=demo.get_basemap()
    fig = Figure(size = (1200, 800), backgroundcolor = :grey80)
    ax = Axis(fig[1, 1])
    im=image!(ax,lon[1,1]..lon[end,1],lat[1,1]..lat[1,end],basemap)
    #hidedecorations!(ax)
    fig,ax
end

"""
    plot(𝐼::Individuals)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.

```
using Pkg; Pkg.activate(temp=true)
Pkg.add.(["CSV", "DataFrames", "FileIO", "Colors", "GLMakie"])

using CSV, DataFrames, FileIO, Colors, GLMakie
include("global_ocean_plotting.jl")

fil=joinpath("inputs","GulfStream_21_27.csv")
df=CSV.read(fil,DataFrame)

nt=length(unique(df.t))
#xlims=(-85.0,5.0); ylims=(20.0,67.0)
xlims=(0.0,360.0); ylims=(-80.0,90.0)
cm=:linear_wcmr_100_45_c42_n256
cr=(-1300,00)

fig,tt=PlottingFunctions.plot([],df,colormap=cm,colorrange=cr)
    #xlims=xlims,ylims=ylims,
fig

file_output_mp4=tempname()*".mp4"
PlottingFunctions.record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end
```
"""
function plot(𝐼,🔴;time=0,xlims=(0.0,360.0),ylims=(-80.0,90.0),
    colormap=:linear_wcmr_100_45_c42_n256,colorrange=(-1300,00),
    add_colorbar=false)

    fig,ax=background()

    np=Int(maximum(🔴.ID))
    nt=length(unique(🔴.t))
    ii=1:10000

    tmp1=🔴[np*0 .+ ii,:lon].!==🔴[np*(nt-1) .+ ii,:lon]
    tmp2=🔴[np*0 .+ ii,:lat].!==🔴[np*(nt-1) .+ ii,:lat]
    jj=ii[findall(tmp1.*tmp2)] 

    time==0 ? tt=Observable(nt) : tt=Observable(time)
    tmp1=groupby(🔴, :t)
    lon_t1=tmp1[1][jj,:lon]
    lat_t1=tmp1[1][jj,:lat]
    
    scatter!(ax,lon_t1,lat_t1,markersize=1.0,color=:lightblue)
    for tx in -12:0 
        ttt=@lift(max(1,$tt+tx))
        lon_tt=@lift(lon360.(tmp1[$ttt][jj,:lon]))
        lat_tt=@lift(tmp1[$ttt][jj,:lat])
        d_tt=@lift(max.(tmp1[$ttt][jj,:d],Ref(-1200)))
        scatter!(ax,lon_tt,lat_tt,markersize=4.0,
        color=d_tt,colorrange=colorrange,colormap=colormap)
    end

    limits!(ax,xlims...,ylims...)

    add_colorbar ? Colorbar(fig[1,2],colorrange=colorrange,colormap=colormap) : nothing

    return fig,tt
end

end #module Plotting 

