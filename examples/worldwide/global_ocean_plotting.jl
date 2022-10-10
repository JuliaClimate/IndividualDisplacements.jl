
module PlottingFunctions 

using GLMakie, IndividualDisplacements, DataFrames

"""
    plot(𝐼::Individuals)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.
"""
function plot(𝐼::Individuals,🔴)
    𝐵=𝐼.𝐷.ODL
    xlims=extrema(𝐵.lon)
    ylims=extrema(𝐵.lat)
    fig=Figure()
    ax=Axis(fig[1,1])
    limits!(ax,-180.0,-90.0,20.0,70.0)
    contour!(ax,𝐵.lon,𝐵.lat,permutedims(𝐵.fld),color=:black,levels=0:0.1:4)

    np=maximum(𝐼.🔴.ID)
    nt=Int(round(size(𝐼.🔴,1)/np))
    if "θ" in names(𝐼.🔴)
        ii=findall((!isnan).(𝐼.🔴[np*0 .+ collect(1:10000),:θ]))
    else
        ii=1:10000
    end
    tmp1=𝐼.🔴[np*0 .+ ii,:lon].!==𝐼.🔴[np*(nt-1) .+ ii,:lon]
    tmp2=𝐼.🔴[np*0 .+ ii,:lat].!==𝐼.🔴[np*(nt-1) .+ ii,:lat]
    jj=ii[findall(tmp1.*tmp2)]

    tt=Observable(nt)
    tmp1=groupby(🔴, :t)
    lon_t1=tmp1[1][jj,:lon]
    lat_t1=tmp1[1][jj,:lat]
    lon_tt=@lift(tmp1[$tt][jj,:lon])
    lat_tt=@lift(tmp1[$tt][jj,:lat])
    
    scatter!(ax,lon_t1,lat_t1,markersize=2.0,color=:lightblue)
    scatter!(ax,lon_tt,lat_tt,markersize=4.0,color=:red)

    return fig,tt
end

end #module Plotting 
