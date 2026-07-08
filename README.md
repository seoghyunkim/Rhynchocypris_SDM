## Description of the data and file structure ##


## The dataset includes 9 files: 3 R code files for each Rhynchocypris species, and 6 GeoPackage files.

1. Files
The GeoPackage files consist of stream network file (stream_network_Rhynchocypris.gpkg), 2 basin boundary files (WKMBBSN.gpkg and WKMSBSN.gpkg), and occurrence segment files for the three species (R.oxycephalus_occ.gpkg, R.steindachneri_occ.gpkg, and R.kumgangensis_occ.gpkg).

2. Code/software
Maximum entropy algorithm (MaxEnt; version 3.4.3), implemented in R (version 4.4.1) with R dismo package (version 1.3.16).
R ENMeval package (version 2.0.5.2) was used to optimize model performance  and tested combinations of feature classes (L, LQ, H, LQH, LQHP, and LQHPT) and regularization multipliers (1–5).

3. Access information
Occurrence data were compiled from the National Ecosystem Survey datasets provided by the National Institute of Ecology and from the National Aquatic Ecological Monitoring Program provided through the Water Environment Information System.
