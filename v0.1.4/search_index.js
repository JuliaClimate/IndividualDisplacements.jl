var documenterSearchIndex = {"docs":
[{"location":"#IndividualDisplacements.jl-1","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.jl","text":"","category":"section"},{"location":"#","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.jl","text":"IndividualDisplacements.jl computes elementary point displacements over a gridded model domain. It can also read / write them from / to file. A typical application is the simulation and analysis of materials drifting or moving over the Global Ocean (e.g. plastics or planktons) or Atmosphere (e.g. dust or chemicals). Inter-operability with popular climate model grids and MeshArrays.jl is an important prospect. IndividualDisplacements.jl was initially designed in relation to MITgcm, ECCOv4 (Forget et al. 2015), and CBIOMES model simulations.","category":"page"},{"location":"#API-Guide-1","page":"IndividualDisplacements.jl","title":"API Guide","text":"","category":"section"},{"location":"#","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.jl","text":"","category":"page"},{"location":"#","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.jl","text":"Modules = [IndividualDisplacements]\nOrder   = [:type,:function]","category":"page"},{"location":"#IndividualDisplacements.ReadDisplacements-Tuple{String,DataType}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.ReadDisplacements","text":"ReadDisplacements(dirIn::String,prec::DataType)\n\nRead displacements from MITgcm output file using MeshArrays and return as a DataFrame.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.ReadGriddedFields-Tuple{}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.ReadGriddedFields","text":"ReadGriddedFields()\n\nRead gridded variables from file using MeshArrays and return result in uvetc Dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.VelComp!-Tuple{Array{Float64,1},Array{Float64,1},Dict,Any}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.VelComp!","text":"VelComp!(du,u,p::Dict,tim)\n\nInterpolate velocity from gridded fields (after exchange on u0,v0) and return position increment du (i.e. x,y,fIndex).\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.VelComp-Tuple{Array{Float64,1},Array{Float64,1},Dict,Any}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.VelComp","text":"VelComp(du,u,p::Dict,tim)\n\nInterpolate velocity from gridded fields and return position increment du\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.VelCopy-Tuple{Any,Any,DataFrames.DataFrame,Any}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.VelCopy","text":"VelCopy(du,u,p::DataFrame,t)\n\nInterpolate velocity from MITgcm float_trajectories output and return position increment du.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.NeighborTileIndices_cs-Tuple{Dict}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.NeighborTileIndices_cs","text":"NeighborTileIndices_cs(grid::Dict)\n\nDerive list of neighboring tile indices for a cs or llc grid + functions that convert indices from one tile to another. Returns a Dict to merge later.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.NeighborTileIndices_dpdo-Tuple{Int64,Int64}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.NeighborTileIndices_dpdo","text":"NeighborTileIndices_dpdo(ni::Int,nj::Int)\n\nList of W, E, S, N neighbor tile IDs in the case of a doubly periodic domain with ni x nj tiles.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.RelocationFunctions_cs-Tuple{MeshArrays.gcmarray}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.RelocationFunctions_cs","text":"RelocationFunctions_cs(xmpl)\n\nDefine matrix of functions to convert indices across neighboring tiles\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.RelocationFunctions_cs_check-Tuple{MeshArrays.gcmarray,Array{Function,2},Int64}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.RelocationFunctions_cs_check","text":"RelocationFunctions_cs_check(xmpl,RF,trgt)\n\nVisualize that RelocationFunctions_cs behaves as expected\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.UpdateLocation_cs!-Tuple{Array{Float64,1},Dict}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.UpdateLocation_cs!","text":"UpdateLocation_cs!\n\nUpdate location (x,y,fIndex) when out of domain. Note: initially, this only works for the dpdo grid type provided by MeshArrays.jl.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.UpdateLocation_dpdo!-Tuple{Array{Float64,1},MeshArrays.gcmgrid}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.UpdateLocation_dpdo!","text":"UpdateLocation_dpdo!\n\nUpdate location (x,y,fIndex) when out of domain. Note: initially, this only works for the dpdo grid type provided by MeshArrays.jl.\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.ex_1-Tuple{}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.ex_1","text":"ex_1()\n\nGlobal ocean case – just reading from file for now.\n\ndf=IndividualDisplacements.ex_1()\n\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"plot_pyplot.jl\"))\nPyPlot.figure(); PlotMapProj(df,300); gcf()\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.ex_2-Tuple{}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.ex_2","text":"ex_2()\n\nReproducing MITgcm/verification/flt_example/ case. This is based on an extended and modified configuration of the standard MITgcm test case.\n\n(df,ref,sol)=IndividualDisplacements.ex_2();\n\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"plot_pyplot.jl\"))\nPyPlot.figure(); PlotBasic(df,300,100000.0); gcf()\n\nusing Plots\nPlots.plot(sol[1,:],sol[2,:],linewidth=5,lc=:black, title=\"One Trajectory Example\",\nxaxis=\"x\",yaxis=\"y\",label=\"Julia Solution\") # legend=false\npl=Plots.plot!(ref[1,:],ref[2,:],lw=3,ls=:dash,lc=:red,label=\"MITgcm Solution\")\n\n\n\n\n\n","category":"method"},{"location":"#IndividualDisplacements.myread-Tuple{String,MeshArrays.gcmarray}","page":"IndividualDisplacements.jl","title":"IndividualDisplacements.myread","text":"myread()\n\nRead a gridded variable from 2x2 tile files. This is used in ReadGriddedFields() with flt_example/\n\n\n\n\n\n","category":"method"}]
}