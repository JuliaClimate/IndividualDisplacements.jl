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

`IndividualDisplacements.jl` [@IndividualDisplacementsFeb2021] is focused on computation and analysis of individual displacements across the global climate system. _Individuals_ is used here as a generic term to represent points, particles, parcels, materials, etc that tend to be carried around the Earth by geophysical fluids (via oceanic currents or atmospheric flows for example). Their _displacements_, or trajectories, are computed by deriving pointwise velocities from flow fields provided by the user and integrating over time (e.g., \autoref{fig:atm} and \autoref{fig:ocn}).

This Julia package is notably aimed at the analysis of numerical models that simulate Oceanic and / or Atmospheric transport processes on staggered C-grids -- those most often used in both regional and global modeling. Inter-operability with global climate model grids like the cube-sphere and lat-lon-cap grid [e.g., @gmd-8-3071-2015] is a key feature enabled by `MeshArrays.jl` [[JuliaCon 2018 (video link)](https://youtu.be/RDxAy_zSUvg), @MeshArraysAugust2020]. The chosen approach also readily supports simpler grids (cartesian, spherical, curvilinear grids) often used in process-oriented models, regional models, or satellite data products.

`IndividualDisplacements.jl` is thus readily suited to exploit climate model output and other gridded data sets in research projects that involve, for example, tracking plankton communities, heat storage, or plastic garbage patches within the Ocean; or dust, water, or chemical compounds within the Atmosphere. To achieve generality and interpretability, `IndividualDisplacements.jl` first defines two data structures (`FlowFields` and `Individuals`) that allow for simple, flexible user specifications of flow fields, initial positions, etc. It then adds a high-level API such that integrating trajectories for individuals $\mathit{I}$ amounts to a single function call, $\int!\,(\,\mathit{I}\,)$.

Internally, the package currently employs `OrdinaryDiffEq.jl` to integrate particle motions over time (within the $\int!\,(\,\mathit{I}\,)$ function call), `NetCDF.jl` and `CSV.jl` for I/O, and `DataFrames.jl` for diagnostics computed along particle trajectories (within $\int!\,(\,\mathit{I}\,)$ or afterwards). The initial test suite [see @IndividualDisplacementsFeb2021] is based not only on idealized flow fields and toy-models, but also on data-constrained Ocean simulations from the OCCA/ECCO projects [@Forget2010, @gmd-8-3071-2015, @Forget2018setup], which are also used in marine ecosystem simulations in CBIOMES [@cbiomes2019]. The OCCA/ECCO gridded ocean circulation estimates are retrieved from a permanent dataverse archive [@OCCAdataverse,@Forget2016dataverse] automatically via the Julia Artifacts system.

`IndividualDisplacements.jl` is also intended to facilitate research involving model-data comparison, data assimilation, or machine learning by providing basic interfaces to other displacement data sets. To start the package provides some support for ingesting (1) data collected in the field by the Global Ocean Drifter Program and Argo array of drifting buoys, and (2) trajectories simulated internally, _online_, by the MIT general circulation model.

The examples folder, which is unit tested upon building the documentation hosted online for each code revision, already covers several common grid cases, two-dimensional and three-dimensional flows, steady-state and time-variable flows. It demonstrates interpolation and diagnostics methods, three plotting libraries, and two notebook systems. Minor parameter changes in these examples should suffice to start applying `IndividualDisplacements.jl` to any other MITgcm simulation of e.g. atmospheric or oceanic turbulence [@MITgcm2020].

![Atmosphere simulation example. Hundred thousand particles are displaced by time-varying 2D flow fields provided by MITgcm on a cube-sphere grid [@MITgcm2020] using `MITgcmTools.jl`. Particle colors show velocity while contours show temperature.\label{fig:atm}](simulated_atm_flow04.png){ width=80% }

![Ocean simulation example. A particle is displaced (left panels) by time-varying 3D flow fields from a global ocean state estimate [@gmd-8-3071-2015] using `OceanStateEstimation.jl`. Vertical profiles of temperature and salinity are recorded along the way (right panels).\label{fig:ocn}](simulated_lagr_float03.png){ width=80% }

# Statement of need

Lagrangian simulation and analysis frameworks such as `IndividualDisplacements.jl` have been widely used across scientific domains involved with geophysical fluids for decades; see [@FLEXPART, @VanSebille2018] for recent reviews in oceanography or atmospheric sciences. `IndividualDisplacements.jl` is most directly related to the `MITgcm/flt` Fortran package which the author also recently extended [@MITgcm2020]. This new Julia package was in fact partly motivated by a need to provide an easier and simpler alternative to using `MITgcm/flt` in its offline mode in order to handle global grids not readily supported by other packages [@gmd-8-3071-2015, @FltGlobalOceanWorkflow2020, @Rousselet2020].

`IndividualDisplacements.jl` notably aims to provide a bridge between the vast community of domain experts, who are typically used to older languages like Fortran, C, Matlab, or Python, and the rapidly growing Julia community and package ecosystem. Julia's native parallelism, GPU support, custom array types, differential programming tools, and plotting libraries indeed all offer great new opportunities for exploiting state of the art HPC climate model simulations and the large data sets that they routinely generate.

The development of `IndividualDisplacements.jl` was motivated not just by our scientific needs and ongoing research projects (see below) but also by the need for more effective tools for climate education, citizen science, and advocacy. In the documentation and examples folder, this is highlighted by the use of unicode in the API, `Pluto.jl` notebooks, and `Makie.jl` animations, which showcase Julia's expressivity and reactivity. 

`IndividualDisplacements.jl` readily supports all common `MITgcm` grids and configurations via `MeshArrays.jl` and `MITgcmTools.jl` (incl. Ocean,  Atmosphere, sea-ice, bio-geo-chemistry, and ecology). This, by itself, yields a large pool of potential scientific applications and expected users; including via the various ongoing research projects that rely on adjoint-optimized, data-constrained solutions from the OCCA/ECCO projects [@Forget2010, @gmd-8-3071-2015]. 

At this stage, `IndividualDisplacements.jl` is considered production-ready, such that we could immediately start transferring several ongoing research collaborations that use MITgcm and ECCO [@Rousselet2020, @Forget2019] from Fortran to Julia. Extension to other models may start with MERRA2 / MITgcm coupled model runs [@Strobach2020] and near-term applications could also include data assimilation to better simulate plastic garbage patches [@gmd-2020-385]. 

The documentation of `IndividualDisplacements.jl` is also now considered sufficient for welcoming additional contributors. Integration with Julia packages developed in `JuliaClimate`, `JuliaOcean`, `JuliaGeo`, `JuliaStats`, `JuliaDynamics`, and other relevant organizations is thus expected to intensify moving forward.

# Acknowledgements

We acknowledge contributions from the open source community at large, the paper reviewers and journal editor, as well as developers of `Julia` and its packages upon which `IndividualDisplacements.jl` crucially depends. 

Funding that supported this work was provided by the Simons Collaboration on Computational Biogeochemical Modeling of Marine Ecosystems (CBIOMES) (grant no. 549931) and the National Aeronautic and Space Administration (grant NASA 19-PO19-0001 from NNH19ZDA001N-PO, grant NASA 80NSSC17K0561 from NNH16ZDA001N-IDS, and grant 19-MAP19-0018 from NNH19ZDA001N-MAP).

# References

