#' Data from observed clones and their properties after electroporation at E12 and sacrifice at P30
#'
#' @format A data frame with 285 rows (clones) and 14 variables:
#' \describe{
#'   \item{Animal_n}{Identifier of the animal the clone was observed in}
#'   \item{clone_ID}{Identifier of the clone, combination of the different colors observed after recombination}
#'   \item{time_of_electroporation_}{Development time when the electroporation was performed}
#'   \item{time_of_analysis}{Time when the sacrifice and analysis were performed}
#'   \item{clone_type}{Type of the clone, either homogeneous (HomC = clone contains cells of the same type), or heterogeneous (HetC = clone contains cells of more than one type)}
#'   \item{clone_subtype}{Subtype of the clone based on the astrocyte types present}
#'   \item{M__L_dispersion_.µm.}{Mediolateral dispersion of the clone, in µm}
#'   \item{clone_size}{Number of cells belonging to the clone (cells expressing the same color combination)}
#'   \item{number_of_BG}{Number of Bergmann Glia cells}
#'   \item{number_of_GLA}{Number of Granular Cell layer Astrocytes}
#'   \item{number_of_WMA}{Number of White Matter Astrocytes}
#'   \item{number_of_CNA}{Number of Cerebellar Nuclei Astrocytes}
#'   \item{number_of_VZA}{Number of Ventricular Zone Astrocyte}
#'   \item{localization}{Comments about the localization of the clone in the cerebellum}
#' }
#' @source Cerrato et al 2018
"Clones_E12_P30"

#' Data from observed clones and their properties after electroporation at E14 and sacrifice at P30
#'
#' @format A data frame with 102 rows (clones) and 14 variables:
#' \describe{
#'   \item{animal_n}{Identifier of the animal the clone was observed in}
#'   \item{clone_ID}{Identifier of the clone, combination of the different colors observed after recombination}
#'   \item{time_of_electroporation_}{Development time when the electroporation was performed}
#'   \item{time_of_analysis}{Time when the sacrifice and analysis were performed}
#'   \item{clone_type}{Type of the clone, either homogeneous (HomC = clone contains cells of the same type), or heterogeneous (HetC = clone contains cells of more than one type)}
#'   \item{clone_subtype}{Subtype of the clone based on the astrocyte types present}
#'   \item{M__L_dispersion_.µm.}{Mediolateral dispersion of the clone, in µm}
#'   \item{clone_size}{Number of cells belonging to the clone (cells expressing the same color combination)}
#'   \item{number_of_BG}{Number of Bergmann Glia cells}
#'   \item{number_of_GLA}{Number of Granular Cell layer Astrocytes}
#'   \item{number_of_WMA}{Number of White Matter Astrocytes}
#'   \item{number_of_CNA}{Number of Cerebellar Nuclei Astrocytes}
#'   \item{number_of_VZA}{Number of Ventricular Zone Astrocyte}
#'   \item{localization}{Comments about the localization of the clone in the cerebellum}
#' }
#' @source Cerrato et al 2018
"Clones_E14_P30"

#' Data from observed clones and their properties after electroporation at E14 and sacrifice at P0
#'
#' @format A data frame with 124 rows (clones) and 14 variables:
#' \describe{
#'   \item{animal_n}{Identifier of the animal the clone was observed in}
#'   \item{clone_ID}{Identifier of the clone, combination of the different colors observed after recombination}
#'   \item{time_of_electroporation_}{Development time when the electroporation was performed}
#'   \item{time_of_analysis}{Time when the sacrifice and analysis were performed}
#'   \item{clone_type}{Type of the clone, either homogeneous (HomC = clone contains cells of the same type), or heterogeneous (HetC = clone contains cells of more than one type)}
#'   \item{clone_subtype}{Subtype of the clone either Cortical (all astrocytes are located in the cortex), WM (all astrocytes in the white matter), Cortical_WMA (astrocytes both in the cortex and in the white matter)}
#'   \item{clone_subtype.1}{Same definition as clone_subtype but here the detail of the cortical astrocyte is given (BG progenitor or GLA)}
#'   \item{M__L_dispersion_.µm.}{Mediolateral dispersion of the clone, in µm}
#'   \item{clone_size}{Number of cells belonging to the clone (cells expressing the same color combination)}
#'   \item{number_of_BG}{Number of Bergmann Glia cells}
#'   \item{number_of_GLA}{Number of Granular Cell layer Astrocytes}
#'   \item{number_of_WMA}{Number of White Matter Astrocytes}
#'   \item{number_of_VZA}{Number of Ventricular Zone Astrocyte}
#'   \item{localization}{Comments about the localization of the clone in the cerebellum}
#' }
#' @source Cerrato et al 2018
"Clones_E14_P0"

#' Data from observed clones and their properties after electroporation at E12 and sacrifice at P0
#'
#' @format A data frame with 63 rows (clones) and 14 variables:
#' \describe{
#'   \item{animal_n}{Identifier of the animal the clone was observed in}
#'   \item{clone_ID}{Identifier of the clone, combination of the different colors observed after recombination}
#'   \item{time_of_electroporation_}{Development time when the electroporation was performed}
#'   \item{time_of_analysis}{Time when the sacrifice and analysis were performed}
#'   \item{clone_type}{Type of the clone, either homogeneous (HomC = clone contains cells of the same type), or heterogeneous (HetC = clone contains cells of more than one type)}
#'   \item{clone_subtype}{Subtype of the clone either Cortical (all astrocytes are located in the cortex), WM (all astrocytes in the white matter), Cortical_WMA (astrocytes both in the cortex and in the white matter)}
#'   \item{clone_subtype.1}{Same definition as clone_subtype but here the detail of the cortical astrocyte is given (BG progenitor or GLA)}
#'   \item{M__L_dispersion_.µm.}{Mediolateral dispersion of the clone, in µm}
#'   \item{clone_size}{Number of cells belonging to the clone (cells expressing the same color combination)}
#'   \item{number_of_BG}{Number of Bergmann Glia cells}
#'   \item{number_of_GLA}{Number of Granular Cell layer Astrocytes}
#'   \item{number_of_WMA}{Number of White Matter Astrocytes}
#'   \item{number_of_VZA}{Number of Ventricular Zone Astrocyte}
#'   \item{folium}{The folium where the clone is located}
#' }
#' @source Cerrato et al 2018
"Clones_E12_P0"

#' Example of parameters used to govern multipotent progenitors behavior
#'
#' @format A list of 3 elements:
#' \describe{
#'   \item{probability}{Matrix with columns: probabilities of generating a MP (1st column), or an Astrocyte (2nd column), 
#'   
#'                                   rows :  for each generation in the simulated lineage (number of rows can vary)}
#'   \item{type}{The corresponding names of different possible outcome, in the same order as the columns of the probability matrix (MP, Astro) 
#'               (number of types correponds to the number of columns of the above probability matrix)}
#'   \item{generationInterval}{Vector of length corresponding to the number of generations, length equal to the number of rows of the matrix probability}
#' }
#' @source Cerrato et al 2018
"transition.MP"


#' Example of parameters used to govern postmitotic astrocyte differentiation produced by a multipotent progenitor
#'
#' @format A list of 3 elements:
#' \describe{
#'   \item{probability}{Matrix with columns: probabilities of generating a BG (1st column), a GLA (2nd column), or a WMA (3rd column)
#'   
#'                                   rows :  for each generation in the simulated lineage (number of rows can vary)}
#'   \item{type}{The corresponding names of different possible outcome, in the same order as the columns of the probability matrix (BG, GLA, WMA) 
#'               (number of types correponds to the number of columns of the above probability matrix)}
#'   \item{generationInterval}{Vector of length corresponding to the number of generations, length equal to the number of rows of the matrix probability}
#' }
#' @source Cerrato et al 2018
"transition.Astro"

#' Example of parameters used to govern the fate of postmitotic astrocyte directly targeted by the electroporation (1 cell clone)
#'
#' @format A list of 2 elements:
#' \describe{
#'   \item{probability}{Vector: probabilities of the one cell clone to be a BG (1st column), a GLA (2nd column), or a WMA (3rd column)}
#'   \item{type}{The corresponding names of different possible outcome, in the same order as the columns of the probability matrix (BG, GLA, WMA) 
#'               (number of types correponds to order of the probability vector element)}
#' }
#' @source Cerrato et al 2018
"firstPostMitoticCell"

#' Example of initialization of the cloneTable for the divisionMPGenerationDpd function
#'
#' @format A dataframe with 4 columns containing initially the data for the 1st mother cell:
#' \describe{
#'   \item{motherID}{ID of the mother of the current cell}
#'   \item{cellID}{ID of the current cell}
#'   \item{type}{Type of the current cell}
#'   \item{timepoint}{Generation the current cell is born. Mother cell is initialized at generation 1}
#' }
#' @source Cerrato et al 2018
"cloneTable"

#' Example of initialization of the first mother cell for the divisionMPGenerationDpd function
#'
#' @format A dataframe with 1row and 4 columns containing the data for the 1st mother cell:
#' \describe{
#'   \item{cellID}{ID of the 1st cell of the lineage}
#'   \item{type}{Type of the 1st cell of the lineage}
#'   \item{timepoint}{Set to generation 1}
#' }
#' @source Cerrato et al 2018
"currentCell"