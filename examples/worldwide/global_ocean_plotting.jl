
module PlottingFunctions 

using GLMakie, IndividualDisplacements, DataFrames
using FileIO, Colors

function background()
    dx=0.1
    path_img="/Users/gforget/mywork/data/MITgcm_highres_sample/images/"
    earth_img=load(joinpath(path_img,
             "Blue_Marble_Next_Generation_topography_bathymetry.jpg"))
    earth_img=reverse(permutedims(earth_img),dims=2)
    fig = Figure(resolution = (1200, 800), backgroundcolor = :grey80)
    ax = Axis(fig[1, 1])
    image!(ax,-179.95:dx:179.95,-89.95:dx:89.95,earth_img)
    hidedecorations!(ax)
    fig,ax
end

"""
    plot(𝐼::Individuals)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.
"""
function plot(𝐼,🔴;time=0,xlims=(-180.0,180.0),ylims=(-90.0,90.0))
    if false
        𝐵=𝐼.𝐷.ODL
        fig=Figure()
        ax=Axis(fig[1,1])
        contour!(ax,𝐵.lon,𝐵.lat,permutedims(𝐵.fld),color=:black,levels=0:0.1:4)
    else
        fig,ax=background()
    end

    np=Int(maximum(🔴.ID))
    nt=length(unique(🔴.t))
    if "θ" in names(🔴)
        ii=findall((!isnan).(🔴[np*0 .+ collect(1:10000),:θ]))
    else
        ii=1:10000
    end
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
          color=d_tt,colorrange=(-1500,0),colormap=:plasma)
    end
    #more time steps

    limits!(ax,xlims...,ylims...)

    return fig,tt
end

end #module Plotting 
