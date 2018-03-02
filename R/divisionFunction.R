#' Simulate an astrocyte lineage
#'
#' @param motherCell the information table of the mother cell
#' @param cloneOutput the initialized clone table with the 1st cell information
#' @param transitionMatrix a list specifying the transition matrix of the particular progenitor type:
#'
#'  - \strong{transition.Progenitor} :  list of 3 elements defining the probabilities of cell cycle exit vs amplification of the progenitor, as well as the potential types produced
#'  \describe{
#'   \item{probability}{Matrix with columns: probabilities of generating a MP (1st column), or an Astrocyte (2nd column), 
#'   
#'                                   rows :  for each generation in the simulated lineage (number of rows can vary)}
#'   \item{type}{The corresponding names of different possible outcome, in the same order as the columns of the probability matrix (MP, Astro) 
#'               (number of types correponds to the number of columns of the above probability matrix)}
#'   \item{generationInterval}{Vector of length corresponding to the number of generations, length equal to the number of rows of the matrix probability}
#' }
#'
#'  - \strong{transition.Astro} : probabilities of generating the different astrocyte types after cell cycle exit
#' \describe{
#'   \item{probability}{Matrix with columns: probabilities of generating a BG (1st column), a GLA (2nd column), or a WMA (3rd column)
#'   
#'                                   rows :  for each generation in the simulated lineage (number of rows can vary)}
#'   \item{type}{The corresponding names of different possible outcome, in the same order as the columns of the probability matrix (BG, GLA, WMA) 
#'               (number of types correponds to the number of columns of the above probability matrix)}
#'   \item{generationInterval}{Vector of length corresponding to the number of generations, length equal to the number of rows of the matrix probability}
#' }
#'
#'  - \strong{firstPostMitoticCell} : for the one-cell clones, proportions of the observed astrocyte types
#'
#' @param maxCount number of trials after which the function will stop and return an error message if the lineage never stops (stochasticity issue, the lineage stops when all 
#'                 cells differentiate. If the probability of cell cycle exit is too low, some cell might keep on proliferating indefinitely.)
#'                 
#' @param parameterInterpolation boolean indicating whether interpolation was performed on the input data (always TRUE in the published simulations)
#'
#' @return cloneTable: dataframe containing the information about the simulated lineage
#' \describe{
#'   \item{motherID}{ID of the mother of the current cell}
#'   \item{cellID}{ID of the current cell}
#'   \item{type}{Type of the current cell}
#'   \item{timepoint}{Generation the current cell is born. Mother cell is initialized at generation 1}
#' }
#'
#' @seealso \code{\link{currentCell}} \code{\link{cloneTable}} \code{\link{transition.MP}} \code{\link{transition.Astro}} \code{\link{firstPostMitoticCell}}
#'
#' @examples
#' # Simulation of a lineage with a multipotent progenitor as initial cell.
#' cloneTable  <- divisionMPGenerationDpd(motherCell = currentCell, 
#'                                        cloneOutput = cloneTable, 
#'                                        transitionMatrix = list(transition.MP, transition.Astro, firstPostMitoticCell),
#'                                        maxCount = 100000)
#'
#' @export divisionMPGenerationDpd
divisionMPGenerationDpd <- function(motherCell = currentCell, cloneOutput = cloneTable,
                                    transitionMatrix = list(transition.Progenitor, transition.Astro, firstPostMitoticCell),
                                    maxCount = 100000,
                                    parameterInterpolation = T){
  # transition.MP : list of probability of becoming MP or Astrocyte, corresponding cell types, generationInterval : for generation dependent probabilities
  # example : transition.Astro <-  list(probability = matrix(c(pBG, pGLA,pWMA), nrow = length(pBG), ncol = 3), type = c("BG", "GLA","WMA"), generationInterval = generationInterval)

  # checks the consistency of the input

  if(parameterInterpolation == F){
    if(dim(as.matrix(transitionMatrix[[2]]$probability))[1] != (length(transitionMatrix[[2]]$generationInterval)-1)
       | dim(as.matrix(transitionMatrix[[2]]$probability))[2] != (length(transitionMatrix[[2]]$type))){
      stop('The dimention of the Astocyte transitions are inconsistent, check the probability, type and generationInterval vectors')
    }
    if(dim(as.matrix(transitionMatrix[[1]]$probability))[1] != (length(transitionMatrix[[1]]$generationInterval)-1)
       | dim(as.matrix(transitionMatrix[[1]]$probability))[2] != (length(transitionMatrix[[1]]$type))){
      stop('The dimention of the MP transitions are inconsistent, check the probability, type and generationInterval vectors')
    }
  } else {
    if(dim(as.matrix(transitionMatrix[[2]]$probability))[1] != (length(transitionMatrix[[2]]$generationInterval))
       | dim(as.matrix(transitionMatrix[[2]]$probability))[2] != (length(transitionMatrix[[2]]$type))){
      stop('The dimention of the Astocyte transitions are inconsistent, check the probability, type and generationInterval vectors')
    }
    if(dim(as.matrix(transitionMatrix[[1]]$probability))[1] != (length(transitionMatrix[[1]]$generationInterval))
       | dim(as.matrix(transitionMatrix[[1]]$probability))[2] != (length(transitionMatrix[[1]]$type))){
      stop('The dimention of the MP transitions are inconsistent, check the probability, type and generationInterval vectors')
    }
    if(length(transitionMatrix[[3]]$probabilities != transitionMatrix[[3]]$type)){
      stop('The dimention of the firstPostMitoticCell probabilitites are inconsistent, check the type and probabilities vectors')
    }
  }
  
  

  count <- 1
  #numbermother <- dim(motherCell)[1]
  while(dim(motherCell)[1]>0){
    if(count > maxCount){
      cloneTable <- rbind(data.frame(motherID = NA, cellID = NA, type = NA, timepoint = NA))
      return(cloneTable)
      stop('Did not converge')
    }

    motherType <- motherCell$type[1];
    motherTimepoint <- motherCell$timepoint[1]; # print(paste("motherType", motherType,"mother cell:", motherCell, "timePoint:",motherTimepoint))
    currentMotherID <- motherCell$cellID[1]

    if(motherType %in% c("WMA", "BG", "Gla")){
      motherCell <- motherCell[-1,] # remove the cell from the list of mother (will not divide)
    }

    if(motherType %in% "Astro"){ # can happen the intial cell in already differentiating
      transition <- transitionMatrix[[3]]
      cellDice <- runif(1)
      cellType <- transition$type[min(which(cellDice < cumsum(transition$probability)))]
      cloneTable$type[count] <- cellType
      motherCell <- motherCell[-1,]
    } else{

      daughterDice <- runif(2)

      transition <- transitionMatrix[[1]]
      progType <- transition$type[1]
      motherGenerationInterval <- .bincode(motherTimepoint, breaks = transition$generationInterval, right = F)

      D1Type <- transition$type[min(which(daughterDice[1] < cumsum(transition$probability[motherGenerationInterval,])))]
      D1ID <- max(cloneTable$cellID) + 1
      D1Tp <- motherCell$timepoint[1]+1
      D2Type <- transition$type[min(which(daughterDice[2] < cumsum(transition$probability[motherGenerationInterval,])))]
      D2ID <-  max(cloneTable$cellID) + 2
      D2Tp <- motherCell$timepoint[1]+1

      #currentMotherID <- motherCell$cellID[1]
      motherCell <- motherCell[-1,]

      if(D1Type == progType){
        motherCell <- rbind(motherCell, data.frame(cellID = c(D1ID), type = c(D1Type), timepoint = D1Tp))
        motherCell$type <- as.character(motherCell$type)
      }else {
        transition <- transitionMatrix[[2]]
        motherGenerationInterval <- .bincode(motherTimepoint, breaks = transition$generationInterval, right = F)
        D1Dice <- runif(1)
        D1Type <- transition$type[min(which(D1Dice < cumsum(transition$probability[motherGenerationInterval,])))]
      }
      if(D2Type == progType){
        motherCell <- rbind(motherCell, data.frame(cellID = c(D2ID), type = c(D2Type), timepoint = D2Tp))
        motherCell$type <- as.character(motherCell$type)
      } else{
        transition <- transitionMatrix[[2]]
        motherGenerationInterval <- .bincode(motherTimepoint, breaks = transition$generationInterval, right = F)
        D2Dice <- runif(1)
        D2Type <- transition$type[min(which(D2Dice < cumsum(transition$probability[motherGenerationInterval,])))]
      }

      cloneTable <- rbind(cloneTable, data.frame(motherID = currentMotherID,
                                                 cellID = c(D1ID, D2ID), type = c(D1Type, D2Type), timepoint = c(D1Tp, D2Tp)))

    }
    count <- count + 1; #print(count); print(cloneTable$timepoint)

  }
  cloneTable
}



#' Add transparency to a color
#'
#' @param col A color.
#' @param alpha The level of transparency.
#' @return The transparent color.
#' @examples
#' color1 <- add.alpha(col = "red", alpha = 0.6)
#'
#' color2 <- add.alpha(colors()[115], 0.5)
#' x <- rnorm(1000, mean = 1, sd = 0.5)
#' y <- rnorm(1000, mean = 1, sd = 0.5)
#' plot(x,y, col = c(color1, color2), pch = 16)
#'
#' @seealso \code{\link{rgb}}, \code{\link{colors}}
#'
#' @export add.alpha
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

getCloneSubtype <- function(cloneOutcome, cellTypes = c("BG", "GLA","WMA")){
  if(sum(cloneOutcome, na.rm = T) >0){
    cloneSubtype <- paste(cellTypes[cloneOutcome > 0 & is.na(cloneOutcome) == F], collapse = "_")
  } else { cloneSubtype <- NA}
  cloneSubtype
}

getCloneType <- function(cloneOutcome){
  cellTypeNumber <- sum(cloneOutcome>0)
  if(!is.na(cellTypeNumber) & cellTypeNumber < 2){
    cloneType <- "HomC"
  } else { cloneType <- "HetC"}
  cloneType
}

#' Add error bars to a bar plot
#'
#' @param x the barplot x coordinates.
#' @param y the data (height of the bars).
#' @param upper matrix of the different error bar values.
#' @param lower matrix of the different error bar values. By default equals to upper
#' @param length length of the horizontal lines of the error bars
#'
#' @examples
#' v1  <- c(1:4)
#' my.x <- barplot(v1)
#' error <- c(0.2, 1, 0.5, 0.1)
#' error.bar(x = my.x, y = v1, upper = error)
#'
#' @export error.bar
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){

  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))

    stop("vectors must be same length")

  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)

}

#' Add error bars to a bar plot
#'
#' @param x the barplot x coordinates.
#' @param y the data (height of the bars).
#' @param upper matrix of the different error bar values.
#' @param lower matrix of the different error bar values. By default equals to upper
#' @param length length of the horizontal lines of the error bars
#'
#' @examples
#' v1  <- c(1:4)
#' my.x <- barplot(v1)
#' error <- c(0.2, 1, 0.5, 0.1)
#' error.bar(x = my.x, y = v1, upper = error)
#'
#' @export meanSDReplicate
## Function to get the mean and SD of the different replates in the model
meanSDReplicate <- function(simulationResults, variable, replicateNumber = 5){

  if(is.null(simulationResults$Replicate)){
    simulationResults$Replicate <- 1
    for ( i in c(1:replicateNumber)){
      simulationResults$Replicate[(1 + (i-1)*numberClonesPerReplicate) : (numberClonesPerReplicate + (i-1) * numberClonesPerReplicate)] <- i
    }
  }
  replicateNumber <- length(unique(simulationResults$Replicate))

  variableCol <- which(names(simulationResults) %in% variable)
  if(length(variableCol) == 0){stop("The variable chosen does not exist in the table, please check spelling and existence")}

  Replicates <- t(matrix(data = unlist(with(simulationResults, tapply(simulationResults[, variableCol], INDEX = Replicate, FUN = table))), ncol = replicateNumber))
  ReplicatesProp <- prop.table(Replicates, margin = 1)

  meanReplicates<- apply(ReplicatesProp, MARGIN = 2, FUN = mean)
  sdReplicates <- apply(ReplicatesProp, MARGIN = 2, FUN = sd)

  list(mean = meanReplicates, sd = sdReplicates)
}

#' Data from observed clones after electroporation at E14 and sacrifice at P30
#'
#' A dataset containing the information of all clones measured
#'
#' @format A data frame with 102 rows and 14 variables:
#' \describe{
#'   \item{animal_n}{Identifier of the animal the clone was observed in}
#'   \item{clone_ID}{Identifier of the clone, combination of the different colors observed after recombination}
#'   \item{time_of_electroporation_}{Development time when the electroporation was performed}
#'   \item{time_of_analysis}{Time when the sacrifice and analysis were performed}
#'   \item{clone_type}{Type of the clone, either homogeneous (HomC = clone contains cells of the same type), or heterogeneous (HetC = clone contains cells of more than one type)}
#'   ...
#' }
#' @source Cerrato et al 2018
"Clones_E14_P30"
