
module PlottingFunctions 

using GLMakie, IndividualDisplacements, DataFrames

"""
    plot(𝐼::Individuals)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.
"""
function plot(𝐼::Individuals,🔴;time=0,xlims=(-180.0,180.0),ylims=(-90.0,90.0))
    𝐵=𝐼.𝐷.ODL
    fig=Figure()
    ax=Axis(fig[1,1])
#    limits!(ax,-180.0,-90.0,20.0,70.0)
    contour!(ax,𝐵.lon,𝐵.lat,permutedims(𝐵.fld),color=:black,levels=0:0.1:4)

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
    lon_tt=@lift(tmp1[$tt][jj,:lon])
    lat_tt=@lift(tmp1[$tt][jj,:lat])
    
    scatter!(ax,lon_t1,lat_t1,markersize=2.0,color=:lightblue)
    scatter!(ax,lon_tt,lat_tt,markersize=4.0,color=:red)

    limits!(ax,xlims...,ylims...)

    return fig,tt
end

end #module Plotting 
