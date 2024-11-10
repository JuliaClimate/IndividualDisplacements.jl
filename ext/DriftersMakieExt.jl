module DriftersMakieExt

	using Makie, Drifters
	import Drifters: InDiPlot, DataFrame, demo
	import Makie: plot

	function plot(x::InDiPlot)
		if !isempty(x.options)
			o=x.options
			if string(o.plot_type)=="simple_plot1"
				simple_plot1(x.data[:I],x.data[:ϕ])
			elseif string(o.plot_type)=="simple_plot2"
				simple_plot2(x.data[:I])
			elseif string(o.plot_type)=="global_plot1"
				global_plot1(x.data[:I],x.data[:df])
			elseif string(o.plot_type)=="plot_start_end"
				plot_start_end(x.data[:I])
			else
				println("unknown option (b)")	
			end
		else
			println("unknown option (a)")
		end
	end

##

"""
    simple_plot1(I,ϕ)

```
using Drifters, CairoMakie
include("basics/random_flow_field.jl")
x=InDiPlot( data=(I=I,ϕ=ϕ), options=(plot_type=:simple_plot1,) )
plot(x)
```
"""
function simple_plot1(I,ϕ)
	I_t = groupby(I, :t)
	nt=length(I_t)

	time = Observable(nt)
	xp=@lift( I_t[$time].x )
	yp=@lift( I_t[$time].y )

	siz=size(ϕ)
	xx=-0.5 .+ collect(1:siz[1])
	yy=-0.5 .+ collect(1:siz[2])
	ll=collect(-1.0:0.2:1.0)*maximum(ϕ)
	
	fig=Figure(size = (900, 600)); set_theme!(theme_light())
	a = Axis(fig[1, 1],xlabel="x",ylabel="y", title="Positions and Streamfunction")		
	contourf!(a, xx,yy, ϕ, levels=ll, colormap=:grayC)
	scatter!(a,I_t[1].x,I_t[1].y,color=:green2)
	scatter!(a,xp,yp,color=:red)

	fig
end

##

"""
    simple_plot2(I)

```
using Drifters, CairoMakie
include("basics/solid_body_rotation.jl")
x=InDiPlot( data=(I=I,), options=(plot_type=:simple_plot2,) )
plot(x)
```
"""
function simple_plot2(I)
	nt=length(I.🔴.x)
	
	time = Observable(nt)
	xx=@lift( [I.🔴.x[1:$time];fill(NaN,nt-$time)] )
	yy=@lift( [I.🔴.y[1:$time];fill(NaN,nt-$time)] )
	zz=@lift( [I.🔴.z[1:$time];fill(NaN,nt-$time)] )
	
	set_theme!(theme_light())
	f=Figure(size = (900, 600))
	a = Axis3(f[1, 1],xlabel="x",ylabel="y",zlabel="z",
		title="Solid body rotation / Spiral example")		
	lines!(a,xx,yy,zz,linewidth=1.0,color=:black)
	scatter!(a,[I.🔴.x[1]],[I.🔴.y[1]],[I.🔴.z[1]],color=:red)
	scatter!(a,[I.🔴.x[nt]],[I.🔴.y[nt]],[I.🔴.z[nt]],color=:green)

	f
end

##

#using Makie, DataFrames, FileIO, Colors
#using DataDeps, MeshArrays, Drifters

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
    global_plot1(I::Individuals)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.

```
using Drifters, GLMakie
include("worldwide/global_ocean_circulation.jl")

x=InDiPlot( data=(I=𝐼,df=tmp_🔴,), options=(plot_type=:global_plot1,) )
fig,tt=plot(x)
fig

file_output_mp4=tempname()*".mp4"
record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end
```
"""
function global_plot1(I::Individuals,🔴::DataFrame;
	time=0,xlims=(0.0,360.0),ylims=(-80.0,90.0),
    colormap=:linear_wcmr_100_45_c42_n256,
	colorrange=(-1300,00),add_colorbar=false)

    fig,ax=background()

    np=Int(maximum(🔴.ID))
    nt=length(unique(🔴.t))
    ii=1:min(10000,np)

    tmp1=🔴[np*0 .+ ii,:lon].!==🔴[np*(nt-1) .+ ii,:lon]
    tmp2=🔴[np*0 .+ ii,:lat].!==🔴[np*(nt-1) .+ ii,:lat]
    jj=ii[findall(tmp1.*tmp2)] 

    🔴_by_t=groupby(🔴, :t)
    time==0 ? tt=Observable(nt) : tt=Observable(time)    
    for tx in -12:0 
        ttt=@lift(max(1,$tt+tx))
        lon_tt=@lift(lon360.(🔴_by_t[$ttt][jj,:lon]))
        lat_tt=@lift(🔴_by_t[$ttt][jj,:lat])
        d_tt=@lift(max.(🔴_by_t[$ttt][jj,:d],Ref(-1200)))
        scatter!(ax,lon_tt,lat_tt,markersize=4.0,
        color=d_tt,colorrange=colorrange,colormap=colormap)
    end

    lon_t1=🔴_by_t[1][jj,:lon]
    lat_t1=🔴_by_t[1][jj,:lat]
    scatter!(ax,lon_t1,lat_t1,markersize=1.0,color=:lightblue)

    limits!(ax,xlims...,ylims...)

    add_colorbar ? Colorbar(fig[1,2],colorrange=colorrange,colormap=colormap) : nothing

    return fig,tt
end

##

"""
    plot_start_end(I::Individuals)

Plot the initial and final positions as scatter plot in `lon,lat` or `x,y` plane.
"""
function plot_start_end(I::Individuals)
🔴_by_t = Drifters.DataFrames.groupby(I.🔴, :t)
set_theme!(theme_light())
fig=Figure(size = (900, 600))
try
	a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
	scatter!(a,🔴_by_t[1].lon,🔴_by_t[1].lat,color=:green2)
	scatter!(a,🔴_by_t[end].lon,🔴_by_t[end].lat,color=:red)
catch
	a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
	scatter!(a,🔴_by_t[1].x,🔴_by_t[1].y,color=:green2)
	scatter!(a,🔴_by_t[end].x,🔴_by_t[end].y,color=:red)
end
return fig
end

end
