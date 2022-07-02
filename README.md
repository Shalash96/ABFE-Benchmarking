# Benchmarking an Absolute Binding Free Energy Calculation Workflow for Computer-Aided Drug Design
---
###### This repo contains required steps and some Python/Bash scripts using to calculate the Absolute Binding Free Energry of Cyclophilin D with 2 fragments and a merged ligand from the data provided by [Bayer AG](https://github.com/bayer-science-for-a-better-life/abfe-benchmark) using [BioSimSpace](https://biosimspace.org/) framework. 
***This done as my MSc graduation project in computational chemistry from the The University Of Edinburgh.***
#
#
#
## Background
###### Over the past five years Free Energy Perturbation (FEP) methodologies have become established in industry to support calculation of binding energies of ligands, primarily at hit to lead and lead optimisation stages. This is particularly the case with relative binding free energy (RBFE) methodologies (Cournia 2017). [Our group](https://github.com/michellab) is developing the software framework [BioSimSpace](www.biosimspace.org) an interoperable software interface for various FEP simulation engines (SOMD, Gromacs, Amber) (Hedges 2019).
###### There is also growing interest in so-called Absolute Binding Free Energy calculation (ABFE) methodologies for Computer-Aided Drug Design (CADD) (Cournia 2020). These are more generally applicable than RBFE methodologies because they do not require a structurally-related reference ligand in the same binding pose. Our group is currently developing an ABFE protocol that relies on the FEP simulation engine SOMD. We are also working on adding an interface to SOMD to BioSimSpace to facilitate high-throughput ABFE calculations. 

###### ABFE methods are still relatively new and haven’t been used as widely as RBFE methodologies. A topic of current interest for the field is to establish the robustness and applicability of ABFE to problems of relevance to CADD. Our group is involved with a nascent community effort spearheaded by the pharmaceutical company Bayer to benchmark various ABFE workflows on similar systems.  Details are available at [Bayer AG]( https://github.com/bayer-science-for-a-better-life/abfe-benchmark). A few studies benchmarking ABFE with single workflows have recently been published (Khalak 2021, Chen 2022, Feng 2022).
