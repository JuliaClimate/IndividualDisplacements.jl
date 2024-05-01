
module PlottingFunctions 

using GLMakie, DataFrames, FileIO, Colors

function background()
    dx=0.1
    path_img="inputs"
    earth_img=load(joinpath(path_img,
             "Blue_Marble_Next_Generation_topography_bathymetry.jpg"))
    earth_img=reverse(permutedims(earth_img),dims=2)
    fig = Figure(size = (1200, 800), backgroundcolor = :grey80)
    ax = Axis(fig[1, 1])
    image!(ax,-179.95:dx:179.95,-89.95:dx:89.95,earth_img)
    hidedecorations!(ax)
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

nt=length(unique(df.t)); xlims=(-85.0,5.0); ylims=(20.0,67.0)
fig,tt=PlottingFunctions.plot([],df,xlims=xlims,ylims=ylims,
	colormap=:linear_wcmr_100_45_c42_n256,colorrange=(-1300,00), add_colorbar=false)
fig

file_output_mp4=tempname()*".mp4"
PlottingFunctions.record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end
```
"""
function plot(𝐼,🔴;time=0,xlims=(-180.0,180.0),ylims=(-90.0,90.0),
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
        lon_tt=@lift(tmp1[$ttt][jj,:lon])
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

