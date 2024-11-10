
module init

using Drifters, MeshArrays, DataFrames, CSV
import Drifters: randn_lonlat

"""
    init_positions(np ::Int)

Randomly distribute `np` points over the Earth, within `P.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_positions(np ::Int; filename="global_ocean_circulation.csv")
    if filename=="global_ocean_circulation.csv"
        p=dirname(pathof(Drifters))
        fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    else
        fil=filename
    end
    return DataFrame(CSV.File(fil))[1:np,:]
end

"""
    init_global_randn(np ::Int , D::NamedTuple)

Randomly distribute `np` points over the Earth, within `D.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , D::NamedTuple)
    (lon, lat) = randn_lonlat(maximum([2*np 10]))
    (_,_,_,_,f,x,y)=InterpolationFactors(D.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(nearest_to_xy(D.msk,x[m],y[m],f[m]).==1.0)[1:np]
    xyf=permutedims([x[m[n]] y[m[n]] f[m[n]]])
    return DataFrame(x=xyf[1,:],y=xyf[2,:],fid=xyf[3,:])
end

"""
    init_gulf_stream(np ::Int , D::NamedTuple)

Randomly distribute `np` points in the Florida Strait region, within 
`D.msk` region, and return position in grid index space (`i,j,subdomain`).
"""
function init_gulf_stream(np ::Int , D::NamedTuple; zs=0:27)
	lons=[-81,-79]
	lats=[26,28]
	lon=rand(2*np)*diff(lons)[1].+lons[1]
	lat=rand(2*np)*diff(lats)[1].+lats[1]
	
	(_,_,_,_,f,x,y)=Drifters.InterpolationFactors(D.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(Drifters.nearest_to_xy(D.msk,x[m],y[m],f[m]).==1.0)[1:np]
    xyf=permutedims([x[m[n]] y[m[n]] f[m[n]]])

	z=zs[1] .+rand(np)*(zs[end]-zs[1])
    return DataFrame(x=xyf[1,:],y=xyf[2,:],z=z,fid=xyf[3,:])
end


"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see the 
global ocean example.
"""
function initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0), level=1)
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])

   return DataFrame(:x => x, :y => y, :z => fill(level,nf),:fid => fill(1,nf))
end

end
