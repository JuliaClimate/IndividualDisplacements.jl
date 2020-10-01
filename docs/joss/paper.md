---
title: 'IndividualDisplacements.jl: a Julia package to simulate and study particle displacements within the climate system'
tags:
  - Julia
  - MITgcm
  - ocean
  - atmosphere
  - Lagrangian
  - particles
  - parcels
  - materials
  - substances
  - tracking
  - pathways
  - trajectories
  - transit times
  - transports
  - dynamics
  - individuals
authors:
  - name: Gael Forget
    orcid: 0000-0003-0872-7098
    affiliation: "1"
affiliations:
 - name: MIT, EAPS
   index: 1
date: 21 September 2020
bibliography: paper.bib
---

# Summary

`IndividualDisplacements.jl` [@IndividualDisplacementsAugust2020] is focused on computation and analysis of individual displacements across the global climate system. _Individuals_ is used here as a generic term to represent points, particles, parcels, materials, etc that tend to be carried around the Earth by geophysical fluids (via oceanic currents or atmospheric flows for example). 

This package is aimed at the analysis of numerical models that simulate Oceanic and / or Atmospheric transport processes on, typically, global staggered C-grids [e.g. see @gmd-8-3071-2015]. Inter-operability with these climate model grids, based on `MeshArrays.jl` which was introduced at [JuliaCon 2018 (video link)](https://youtu.be/RDxAy_zSUvg) [@MeshArraysAugust2020], is a key feature. The chosen approach covers the comparatively simple grids used in e.g. satellite data products as well.

`IndividualDisplacements.jl` is readily suited to exploit climate model output and other gridded data sets in research projects that involve tracking plankton communities, heat storage, or plastic garbage patches within the Ocean for example; dust, water, or chemical compounds within the Atmosphere, etc. It currently employs `OrdinaryDiffEq.jl` to integrate particle motions over time, `NetCDF.jl` and `CSV.jl` for I/O, and `DataFrames.jl` for diagnostics computed along particle trajectories.

The package is also intended to facilitate research involving, for example, model-data comparison, model intercomparison, data assimilation, state estimation, or machine learning by providing basic interfaces to observed displacement data sets (initially those collected by the Global Ocean Drifter Program and Argo array of drifting buoys) and to simulated displacement data sets generated internally, _online_, by climate models, or via other Lagrangian toolboxes (initially those from the MIT general circulation model).

The initial test suite is based on data-constrained Ocean simulations from the OCCA/ECCO project [@Forget2010],[@gmd-8-3071-2015], [@Forget2018setup], also used in marine ecosystem simulations in CBIOMES [@cbiomes2019], along with idealized geophysical turbulence flow fields [see @IndividualDisplacementsAugust2020]. The OCCA/ECCO gridded ocean circulation estimates are retrieved from a permanent dataverse archive [@OCCAdataverse], [@Forget2016dataverse].

The examples folder, which is unit tested upon building the documentation, already covers several common global grids, cases of two-dimensional and three-dimensional flows, steady-state and time-variable flows, as well as interpolation and diagnostics methods, three plotting libraries, and two notebook systems. Minor parameter changes in these examples should suffice to start applying `IndividualDisplacements.jl` to any other MITgcm simulation of e.g. atmospheric or oceanic turbulence [@MITgcm2020].

# Statement of need 

Lagrangian simulation and analysis frameworks such as `IndividualDisplacements.jl` have been widely used across scientific domains involved with geophysical fluids for decades [see @FLEXPART],[@VanSebille2018 for recent reviews in oceanography or atmospheric sciences]. 

`IndividualDisplacements.jl` provides a bridge between this vast community of domain experts, who are typically used to older languages like Fortran, C, Matlab, or Python, and the rapidly growing Julia community and package ecosystem. Julia's native parallelism, GPU support, custom array types, differential programming tools, and plotting libraries indeed all offer great new opportunities for exploiting state of the art HPC climate model simulations and the large data sets that they routinely generate.

The development of `IndividualDisplacements.jl` was motivated not just by our scientific needs and ongoing research projects (see below) but also by the need for more effective tools for climate education, citizen science, and advocacy. In the documentation and examples folder, this is highlighted by the use of unicode in the API, `Pluto.jl` notebooks, and `Makie.jl` animations, which showcase Julia's expressivity and reactivity. 

`IndividualDisplacements.jl` readily supports all common `MITgcm` grids and configurations via `MeshArrays.jl` (incl. Ocean,  Atmosphere, sea-ice, bio-geo-chemistry, and ecology). This, by itself, yields a large pool of potential scientific applications and expected users -- including via the many ongoing research projects that rely on adjoint-optimized, data-constrained solutions from the OCCA/ECCO project [see @Forget2010][@gmd-8-3071-2015]. 

At this stage, `IndividualDisplacements.jl` is considered production-ready, such that we could immediately start transferring several ongoing research collaborations that use MITgcm and ECCO [e.g. @Rousselet2020], [@Forget2019] from fortran to Julia. Extension to other models may start with MERRA2 / MITgcm coupled model runs later this year [@Strobach2020]. 

The documentation of `IndividualDisplacements.jl` is also now considered sufficient for welcoming additional contributors. Integration with Julia packages developed in `JuliaClimate`, `JuliaOcean`, `JuliaGeo`, `JuliaStats`, `JuliaDynamics`, and other relevant organizations is thus expected to intensify over the coming months.

# Acknowledgements

We acknowledge contributions from the open source community at large, the paper reviewers and journal editor, as well as developers of `Julia` and its packages upon which `IndividualDisplacements.jl` crucially depends. 

Funding that supported this work was provided by the Simons Collaboration on Computational Biogeochemical Modeling of Marine Ecosystems (CBIOMES) (grant no. 549931) and the National Aeronautic and Space Administration (grant NASA 19-PO19-0001, NNH19ZDA001N-PO, and grant NASA 80NSSC17K0561, NNH16ZDA001N-IDS).

# References

