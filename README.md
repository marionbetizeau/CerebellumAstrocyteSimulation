# CerebellumAstrocyteSimulation
Code and package (SimLinTool) to performe astrocyte lineage simulation

This simulation requires R and Rstudio. For downloads see https://cran.r-project.org/ and https://www.rstudio.com/products/rstudio/download/

The simulation functions are in the SimLinTool package , the script SimulationScriptAstrocyteCerebellum.R performs the simulations shown in Cerrato et al.

The datasets are embedded in the SimLinTool package and also available as tab delimited tables in the folder data_tabDelimited

To install the package SimLinTool and run simulations:
1. Install devtools  

          install.pacakges("devtools")
2. Load devtools 

          library(devtools) # NB: this takes a little while
3. Install SimLinTool

          install_github(repo = "marionbetizeau/CerebellumAstrocyteSimulation", auth_token = "593b6d1fc8e79309d70b29ae58b6f7547d47c035")
4. Datasets are: Clones_E12_P0, Clones_E12_P30, Clones_E14_P0, Clones_E14_P30 to get a description of the measured datasets use:  

          ?Clones_E12_P0 
5. Download the script SimulationScriptAstrocyteCerebellum.R, possibility to choose the output folder and change the parameters (see script)
6. Run the entire script using the "source" button top right of Rstudio editor window
