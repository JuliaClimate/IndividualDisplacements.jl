using Test, Documenter
using IndividualDisplacements, Climatology, MeshArrays

Climatology.get_ecco_velocity_if_needed()
Climatology.get_occa_velocity_if_needed()
MeshArrays.GridLoad(MeshArrays.GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
MeshArrays.GridLoad(MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360))
IndividualDisplacements.flt_example_download()

@testset "global" begin
    p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/worldwide/ECCO_FlowFields.jl"))
    𝑃,𝐷=ECCO_FlowFields.global_ocean_circulation()
    df = ECCO_FlowFields.init_from_file(10)
    𝐼=Individuals(𝑃,df.x,df.y,df.f,(𝐷=𝐷,))
    𝑇=(0.0,𝐼.𝑃.𝑇[2])
    ∫!(𝐼,𝑇)

    add_lonlat!(𝐼.🔴,𝐷.XC,𝐷.YC)
    add_lonlat!(𝐼.🔴,𝐷.XC,𝐷.YC,𝑃.update_location!)
    tmp=interp_to_xy(𝐼.🔴,𝐷.YC)
    gcdist(𝐼)

    @test prod(abs.(tmp).<90.0)

    tmp1=randn_lonlat(10)
    tmp2=stproj_inv(stproj(30.0,30.0)...)
    @test prod(isapprox.(tmp2,30.0,atol=1.0))
end

@testset "various" begin
    u,v,w,pos=random_flow_field(format=:Array)
    𝐹=FlowFields(u,u,v,v,[0,1.0])
    𝐼=Individuals(𝐹,pos...)
    ∫!(𝐼)
    
    show(𝐼)
    diff(𝐼)
    size(𝐼)
    𝐽=similar(𝐼)
    @test isa(𝐽,Individuals)

    𝐺=convert_to_FlowFields(u,v,10.0)
    tmp2=nearest_to_xy(𝐺.u0,3.,3.,1.)
    @test isa(tmp2,Array)
    tmp3=nearest_to_xy(𝐹.u0,3.,3.)
    @test isa(tmp3,Array)
end

@testset "doctests" begin
    doctest(IndividualDisplacements; manual = false)
end
