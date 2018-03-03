# author: "Marion Betizeau"
# date: "February 23, 2018"
# Script to perform simulations Figure 9. From Cerrato et al.

# Instructions ------------------------------------------------------------

# 0. Go through the readMe on the GitHub repository to install R, RStudio
# i. Chose the parameter set in the sections: 1, 2, 3 and 4
# ii. Run the entire script click on the Source button in the tool bar (top right of RStudio Editing Window)


# 1. Global Parameters --------------------------------------------------------------

# Folder to save the results
directory <- getwd() 
cat(paste("Results will be stored in", directory))
cat("\n")
cat("To change the output directory use the set setwd, info with ?setwd");cat("\n")

# Loads the required R packages
library(LinSimTool) # Custom made package available on the gitHub repository to install:    install_github(repo = "marionbetizeau/CerebellumAstrocyteSimulation", auth_token = "593b6d1fc8e79309d70b29ae58b6f7547d47c035")
library(igraph)

# Colors associated with BG, GLA, WMA
myColors <- c("green", "yellow", "red")

# Choose the desire age and cerebellar region
age <-    "E12"       # or "E14"
region <- "hemisphere" # or "vermis"  
cat("\n")
cat("\n")
cat(age, region);cat("\n") # print the parameters in the console

# Frequencies of the different progenitor types. 
# MP : multipotent progenitor (For only one multipotent progenitor set freqMP = 1)
# WMA.P : WMA progenitor, produces exclusively WMA
freqMP <-     1
freqWMA.P <-  1 - freqMP  
cat(freqMP*100, "% of MP");cat("\n")

if(freqWMA.P ==0){ # keep track if only 1 or 2 progenitors
  progType <- "MP"
}else{ progType <- "MP_WMAP"}
cat(progType, "progenitor types");cat("\n")

useDataProlif   <- F    # TRUE = using the measured cell cycle exit (from Leto et al, 2011), FALSE = using custom values  (default)
useDataDiff     <- T    # TRUE = using the measured birth dating values, FALSE = using some custom values (see next section to adapt) 
useDataFirst    <- T    # TRUE = using the proportion of single cell clones for the cell that are directly postmitotic (default)

# Number of replicates
replicateNumber <- 5
if(age == "E12") { numberClonesPerReplicate <- 275} 
if(age == "E14") { numberClonesPerReplicate <- 98} 
iterationNumber <- numberClonesPerReplicate * replicateNumber
cat(paste(replicateNumber,"replicates of",numberClonesPerReplicate, "clones"));cat("\n")

# Probabilities that the first cell of a lineage is a progenitor or postmitotic (proportion of one-cell clones in the observed dataset)
if (age == "E12"){
  pInitialCell <- 0.7; # proportions of 1 cell clones in E12-P30 dataset 
  # 0.85 # alternative from Florio et al 2012 growth fraction in VZ at E12
}
if (age == "E14"){
  pInitialCell <- 0.6;  # proportions of 1 cell clones in E14-P30 dataset 
}
cat(pInitialCell*100,"% Progenitors in the starting population");cat("\n")

# 2. Transition Probabilities for Multipotent progenitors (MP) --------------------------------------------------------------

## Here change the parameters for probabilities of generating the different astrocyte types
## Depending if useDataDiff equals T or F, depending on the age and the region go the the right position in the "if statements"
## The generationInterval correspond to the developmental times. 0 is the beginning of the lineage either E12 or E14. 
##      Generation 6 is P0 for E12 starting clones
##      Generation 4 corresponds to P0 for E14 starting clones


if(useDataDiff == T){ # Case using the birth-dating measures
  if(age == "E12"){
    if (region == "vermis"){
      # Birthdating data E12 / vermis 
      generationInterval.Astro <- c(0, 2, 6, 12, 20, 150) # E12, E14, P0, P7, P15+
      #         E12  E14   P0    P7     P15+
      pBG <- c(0.25, 0.33, 0.37, 0.695, 0.83)
      pGLA <- c(0.25, 0.34, 0.4, 0.29, 0.17)
      pWMA <- c(0.5, 0.33, 0.23, 0.015, 0)
    }
    
    if(region == "hemisphere"){
      # Birthdating data E12 / hemisphere 
      generationInterval.Astro <- c(0, 2, 6, 12, 20, 150) # E12, E14, P0, P7, P15+
      #         E12  E14   P0    P7   P15+
      pBG <- c(0.1, 0.1, 0.18, 0.82, 0.87)
      pGLA <- c(0.1, 0.15, 0.36, 0.18, 0.13)
      pWMA <- c(0.8, 0.75, 0.46, 0, 0)
    }
  }
  
  if (age == "E14"){
    if (region == "vermis"){
      # Birthdating data E14 / vermis 
      generationInterval.Astro <- c(0, 4, 8, 12, 150) # E14, P0, P7, P15+
      #         E14  P0    P7    P15+
      pBG <- c(0.52, 0.37, 0.70, 0.83)
      pGLA <- c(0.15, 0.4, 0.29, 0.17)
      pWMA <- c(0.33, 0.23, 0.014, 0) 
      
    }
    if(region == "hemisphere"){
      # Birthdating data E14 / hemisphere 
      generationInterval.Astro <- c(0, 4, 8, 12, 150) # E14, P0, P7, P15+
      #         E14  P0    P7    P15+
      pBG <- c(0.1, 0.22, 0.82, 0.87)
      pGLA <- c(0.15, 0.32, 0.18, 0.13)
      pWMA <- c(0.75, 0.46, 0, 0)
    }
  }
  
  diffParam <- "data"
  
}else{ # Case using custom measures
  
  if (age == "E12"){
    if (region == "vermis"){
      generationInterval.Astro <- c(0, 3, 6, 10, 20, 150) # E12, E14, P0, P7, P15+
      #         E12  E14   P0    P7   P15+
      pBG <- c(0.15, 0.1, 0.22, 0.82, 0.87)
      pGLA <- c(0.15, 0.15, 0.32, 0.18, 0.13)
      pWMA <- c(0.7, 0.75, 0.46, 0, 0)
    }
    if (region == "hemisphere"){
      generationInterval.Astro <- c(0, 2, 6, 10, 20, 150) # E12, E14, P0, P7, P15+
      #         E12  E14   P0    P7   P15+
      pBG <- c(0.05, 0.05, 0.22, 0.82, 0.87)
      pGLA <- c(0.05, 0.05, 0.32, 0.18, 0.13)
      pWMA <- c(0.9, 0.9, 0.46, 0, 0)
    }
  }
  
  if (age == "E14"){
    if (region == "vermis"){
      
      generationInterval.Astro <- c(0, 4, 10, 14, 150) # E14, P0, P7, P15+
      # Observed vermis data E15, P1, P7, P15 (E15 data high variations: take 1/3 of each type)
      #         E14  P0    P7    P15+
      pBG  <- c(0.33, 0.37, 0.70, 0.83)
      pGla <- c(0.34, 0.40, 0.29, 0.17)
      pWMA <- c(0.33, 0.23, 0.014, 0)
    }
    if (region == "hemisphere"){
      generationInterval.Astro <- c(0, 4, 10, 14, 150) 
      # Test values E14 hemisphere
      #         E14  P0    P7    P15+
      pBG <- c(0.1, 0.19, 0.82, 0.87)
      pGLA <- c(0.2,  0.35, 0.18, 0.13)
      pWMA <- c(0.7,  0.46, 0, 0)
    }
  }
  diffParam <- "custom"
}


# 3. Transition Probabilities for WMA progenitors (WMA.P, if it applies) --------------------------------------------------------------

# Production of only WMA
generationInterval.WMAP <- c(0, 15, 150) 
pBG.WMAP <- c(0, 0) 
pGLA.WMAP <- c(0, 0)
pWMA.WMAP <- c(1, 1)



# 4. Cell cycle exit parameters --------------------------------------------------------------

if(useDataProlif == T){ # Case using the measured cell cycle exit based on the results from Leto et al 2011 in the White Matter
  if(age == "E12"){
    # Cell cycle exit based on the results from Leto et al 2011 in the White Matter
    generationInterval.MP <- c(0, 3, 6, 10, 20, 150) # E12, E14, P0, P7, P15+
      #                   E12  E14  P1    P5   P15+
    pamplifMP <- 1 - c(0.55, 0.7, 0.85, 0.5, 1)
    
    # WMA committed progenitor : stop proliferating at around P2
    generationInterval.WMAP <- c(0,6,8,150) # E12, P0, P2+
    #                 E12  P0  P2+
    pamplif.WMAP <- c(0.1, 0.1, 0)
  }
  
  if(age == "E14"){
    # Cell cycle exit based on the results from Leto et al 2011 in the White Matter
    generationInterval.MP <- c(0, 4, 10, 14, 150) # E14, P0, P7, P15+
    #                  E14  P1    P5   P15+
    pamplifMP <- 1 - c(0.7, 0.85, 0.5, 1)
    
    # WMA committed progenitor : stop proliferating at around P2
    generationInterval.WMAP <- c(0,4,6,150) # E14, P0, P2+
    #                 E14  P0 P2+
    pamplif.WMAP <- c(0.2, 0, 0)
  }
  prolifParam <- "data"
  
} else {  # Case using custom measures, default used in the simulation: - constant over cerebellar development for MP, 
          #                                                             - WMA.MP stops proliferating at around P2
  
  if(age == "E12"){
    pamplifMP <- c(0.465,0.465) # Optimal value to reproduce clone size at P30
    generationInterval.MP <- c(0,10,150)
    
    # WMA committed progenitor
    generationInterval.WMAP <- c(0,6,8,150) # E12, P0, P2+
    #                 E12  P0  P2+
    pamplif.WMAP <- c(0.1, 0.1, 0)
  }
  if(age == "E14"){
    pamplifMP <- c(0.44,0.44) # Optimal value to reproduce clone size at P30
    generationInterval.MP <- c(0,10,150)
    
    # WMA committed progenitor
    generationInterval.WMAP <- c(0,4,6,150) # E14, P0, P2+
    #                 E14  P0 P2+
    pamplif.WMAP <- c(0.2, 0, 0)
    
  }
  prolifParam <- "custom"
}



# 5. Load the relevant observed datasets --------------------------------------------------------------

if(age == "E12"){
  observedData <- Clones_E12_P30 # loads the corresponding dataset
}
if(age == "E14"){  
  observedData <- Clones_E14_P30 # loads the corresponding dataset
  observedData$clone_subtype[observedData$clone_subtype == "WMA___CNA"] <- "WMA"
}

# ignore clones with only CA astrocytes
observedData <- observedData[observedData$clone_subtype != "CNA",]
observedData$clone_subtype <- droplevels(observedData$clone_subtype)
observedData$clone_subtype_NoCNA <- apply(observedData[,9:11], 1, FUN = getCloneSubtype)
observedData$clone_type_NoCNA <- apply(observedData[,9:11], 1, FUN = getCloneType)

observedData$clone_size <- observedData$clone_size - observedData$number_of_CNA  # ignores the CNA cells, 

# Computes the ratios of astrpcyte types per clone
totalAstro <- sum(observedData$clone_size)
observedData$GLAoverWMA <- observedData$number_of_GLA / observedData$number_of_WMA
observedData$BGoverWMA <- observedData$number_of_BG / observedData$number_of_WMA
observedData$BGoverGLA <- observedData$number_of_BG / observedData$number_of_GLA



# 6. Parameters for targeted cells that are postmitotic (one cell clones) --------------------------------------------------------------
  # observed proportions of one-cell clones of different types


uniCellclones <- observedData[observedData$clone_size == 1,]
uniCellclones$clone_subtype[uniCellclones$clone_subtype == "WMA_CNA"] <- "WMA"
uniCellclones$clone_subtype <- droplevels(uniCellclones$clone_subtype)
firstp <- prop.table(table(uniCellclones$clone_subtype))

firstPostMitoticCell <- list(probability = matrix(firstp , nrow = 1, ncol = 3), type = c("BG", "GLA","WMA"))

firstPostMitoticCell.WMAP <- list(probability = matrix(c(0,0,1) , nrow = 1, ncol = 3), type = c("BG", "GLA","WMA"))



# 7. Interpolate linearly and plots the input parameters --------------------------------------------------------------

# Interpolate linearly between the different datapoints and plots the input parameters
# For interpolation need to repeat the last value of the parameters
pBG <- c(pBG, pBG[length(pBG)])
pGLA <- c(pGLA, pGLA[length(pGLA)])
pWMA <- c(pWMA, pWMA[length(pWMA)])
pamplifMP <- c(pamplifMP, pamplifMP[length(pamplifMP)])
generations <- seq(min(generationInterval.Astro), max(generationInterval.Astro),1)

generationpBG.WMAP <- rep(pBG.WMAP[1], times = length(generations))
generationpGLA.WMAP <- rep(pGLA.WMAP[1], times = length(generations))
generationpWMA.WMAP <- rep(pWMA.WMAP[1], times = length(generations))
pamplif.WMAP <- c(pamplif.WMAP, pamplif.WMAP[length(pamplif.WMAP)])


if(age == "E12"){
  myXlab <- "Generation, generation 6 ~ P0"
  P0Gen <- 6
}
if(age == "E14"){
  myXlab <- "Generation, generation 4 ~ P0"
  P0Gen <- 4
}
generations <- seq(min(generationInterval.Astro), max(generationInterval.Astro),1)
generationpBG <- rep(NA, length(generations))
generationpGLA <- rep(NA, length(generations))


generationpamplifMP <- rep(NA, length(generations))

if(progType == "MP_WMAP"){
  generationpamplif.WMAP <- rep(NA, length(generations))
  generationpWMA <- rep(NA, length(generations))
  for(i in c(1:(length(pamplif.WMAP)-1))){
    x <-  generationInterval.WMAP[i:(i+1)]
    yamplif.WMAP <- pamplif.WMAP[i:(i+1)]
    my.lmamplif.WMAP <- lm(yamplif.WMAP~x)
    generationpamplif.WMAP[seq(generationInterval.WMAP[i], generationInterval.WMAP[i+1], 1)+1] <- seq(generationInterval.WMAP[i], generationInterval.WMAP[i+1], 1) * my.lmamplif.WMAP$coefficients[2] + my.lmamplif.WMAP$coefficients[1]
  }
}

# Interpolate between transition probabilities, MP and WMA.P
for(i in c(1:(length(pBG)-1))){
  x <-  generationInterval.Astro[i:(i+1)]
  
  yBG <- pBG[i:(i+1)]
  my.lmBG <- lm(yBG~x)
  generationpBG[seq(generationInterval.Astro[i], generationInterval.Astro[i+1], 1)+1] <- seq(generationInterval.Astro[i], generationInterval.Astro[i+1], 1) * my.lmBG$coefficients[2] + my.lmBG$coefficients[1]
  
  yGLA <- pGLA[i:(i+1)]
  my.lmGLA <- lm(yGLA~x)
  generationpGLA[seq(generationInterval.Astro[i], generationInterval.Astro[i+1], 1)+1] <- seq(generationInterval.Astro[i], generationInterval.Astro[i+1], 1) * my.lmGLA$coefficients[2] + my.lmGLA$coefficients[1]
  generationpWMA <- 1 - generationpBG - generationpGLA
  
  if(progType == "MP_WMAP"){
    yWMA <- pWMA[i:(i+1)]
    my.lmWMA <- lm(yWMA~x)
    generationpWMA[seq(generationInterval.WMAP[i], generationInterval.WMAP[i+1], 1)+1] <- seq(generationInterval.WMAP[i], generationInterval.WMAP[i+1], 1) * my.lmWMA$coefficients[2] + my.lmWMA$coefficients[1]
  }
  
}

# Interpolate between amplification probabilities, MP 
for(i in c(1:(length(pamplifMP)-1))){
  x <-  generationInterval.MP[i:(i+1)]
  yamplifMP <- pamplifMP[i:(i+1)]
  my.lmamplifMP <- lm(yamplifMP~x)
  generationpamplifMP[seq(generationInterval.MP[i], generationInterval.MP[i+1], 1)+1] <- seq(generationInterval.MP[i], generationInterval.MP[i+1], 1) * my.lmamplifMP$coefficients[2] + my.lmamplifMP$coefficients[1]
}

# Structure the input parameters
# Input parameters for each cell type (MP, WMA.P and differentiating astrocytes) are lists of
#     - probability: matrix with columns: possible outcomes for this type (ex: MP and astrocyte for an MP progenitor), 
#                                 rows : probabilities of producing the different possible outcomes for this type at each generation 
#     - type: the corresponding names of different possible outcome (number of types correponds to the number of columns of the above probability matrix)
#     - generationInterval : vector of the number of generations, length equal to the number of rows of the matrix probability

transition.Astro <- list(probability = matrix(c(generationpBG, generationpGLA, generationpWMA), nrow = length(generationpBG), ncol = 3),
                         type = c("BG", "GLA","WMA"), generationInterval = generations)
transition.MP <- list(probability = matrix(c(generationpamplifMP, 1 - generationpamplifMP), nrow = length(generationpamplifMP), ncol = 2), 
                      type = c("MP", "Astro"), generationInterval = generations)

transition.Astro.WMAP <- list(probability = matrix(c(generationpBG.WMAP, generationpGLA.WMAP, generationpWMA.WMAP), nrow = length(generationpBG), ncol = 3), 
                              type = c("BG", "GLA","WMA"), generationInterval = generations)

if(progType == "MP_WMAP"){
  generationpAstro.WMAP <- 1 - generationpamplif.WMAP
  transition.WMAP <- list(probability = matrix(c(generationpamplif.WMAP, 1 - generationpamplif.WMAP), nrow = length(generationpamplifMP), ncol = 2), 
                          type = c("WMAP", "Astro"), generationInterval = generations)
}

  



# 8. Simulate the clones and output the lineages as PDF file ------------------------------------------------------------

simulationName <- paste(Sys.time(), progType, age,region,"Diff", diffParam,"Prolif", prolifParam, round(mean(pamplifMP),3),
                        replicateNumber, "replicates_of", numberClonesPerReplicate,"lineages", sep = "_") 

simulationResults <- NULL # initialize the output table
set.seed(12) # to generate the same random number, useful to debug

pdf(file = paste(simulationName,"_lineages.pdf", sep = ""), width = 8, height = 8) # to save the generated lineage

# Simulate the different lineages (iterationNumber times)
for(iter in c(1:iterationNumber)){
  
  # Initialize the type of the first cell
  initCell <- runif(1)
  initType <- c("MP","WMAP","Astro")[min(which(initCell < cumsum(c(pInitialCell * freqMP, pInitialCell * freqWMA.P, 1 - pInitialCell))))] # Possible types of the initial cell
  
  # Initialise the current cell (first cell of the lineage)
  currentCell <- data.frame(cellID = 1, type = initType, timepoint = 1) 
  currentCell$type <- as.character(currentCell$type)
  
  # Initialise the output of the clone
  cloneTable <- data.frame(motherID = NA, cellID = 1, type = initType, timepoint = 1)
  cloneTable$type <- as.character(cloneTable$type)
  
  # Generates clones from the MP and the WMAP if any
  if(initType == "MP" | initType == "Astro"){
    cloneTable  <- divisionMPGenerationDpd(motherCell = currentCell, 
                                           cloneOutput = cloneTable, 
                                           transitionMatrix = list(transition.MP, transition.Astro, firstPostMitoticCell),
                                           maxCount = 10000)
  }else if(initType == "WMAP"){
    cloneTable  <- divisionMPGenerationDpd(motherCell = currentCell, 
                                           cloneOutput = cloneTable, 
                                           transitionMatrix = list(transition.WMAP, transition.Astro.WMAP, firstPostMitoticCell.WMAP),
                                           maxCount = 10000)
  }
  
  # Deals with the replicates:
  
  if (iter ==1){
    storeCloneTable <- cbind(cloneID = iter, initType = initType, cloneTable)
  }else{storeCloneTable <- rbind(storeCloneTable, cbind(cloneID = iter, initType = initType, cloneTable))}
  
  
  # Check the lineage and plot them
  if(nrow(cloneTable) == 1 & is.na(cloneTable$type[1])){
    simulationResults <- rbind(simulationResults, NA)
  }
  else {
    is.leave <- !cloneTable$cellID %in% cloneTable$motherID
    leavesType <- factor(cloneTable$type[is.leave], levels = c("BG", "GLA","WMA"))
    simulationResults <- rbind(simulationResults, table(leavesType))
  }
  
  colorList <- c("gray","coral", myColors)
  
  # Plots the lineage if at least a division
  types <- c("MP","WMAP","BG", "GLA", "WMA")
  
  if(nrow(cloneTable) > 2){ # If at least a division in the lineage
    # Visualize the outcome using the igraph library
    clone.net <- graph.data.frame(cloneTable[2:nrow(cloneTable), ], directed=T)
    myLayout <- layout.reingold.tilford(clone.net)
    
    verticesID <- match(attributes(V(clone.net))$names, cloneTable$cellID)
    
    V(clone.net)$color <- colorList[match(cloneTable$type[verticesID], types)]
    
    # Plot the lineage
    plot(clone.net, layout = myLayout, vertex.size = 5, main = paste(sum(is.leave), "leaves"))
    legend("topleft", legend = types, pch=21,
           pt.bg=colorList, pt.cex=1.5, bty = 'n', cex = 0.8)
  } else{ # If only one cell in the lineage
    colorNumber <- match(cloneTable$type, types)
    plot(0,0, pch = 19 , col = colorList[colorNumber], cex = 2.5, main = "1 leave", xlab = "", ylab = "")
  }
}

dev.off() # closes the pdf device to save the 

# Analysis of the simulated lineages
simulationResults <- as.data.frame(simulationResults)
simulationResults$CloneSize <- rowSums(simulationResults)

# Create clone subtype column
simulationResults$cloneSubtype <- apply(simulationResults[,1:3], 1, FUN = getCloneSubtype)
simulationResults$cloneSubtype <- factor(simulationResults$cloneSubtype, levels = c("BG","BG_GLA","BG_GLA_WMA","BG_WMA","GLA","GLA_WMA","WMA"))
simulationResults$cloneType <- apply(simulationResults[,1:3], 1, FUN = getCloneType)
simulationResults$cloneType <- factor(simulationResults$cloneType, levels = c("HetC", "HomC"))

# Gets the ratio of astrocyte types per lineage
simulationResults$GLAoverWMA <- simulationResults$GLA / simulationResults$WMA
simulationResults$GLAoverWMA[!is.finite(simulationResults$GLAoverWMA)] <- NA

simulationResults$BGoverWMA <- simulationResults$BG / simulationResults$WMA
simulationResults$BGoverWMA[!is.finite(simulationResults$BGoverWMA)] <- NA

simulationResults$BGoverGLA <- simulationResults$BG / simulationResults$GLA
simulationResults$BGoverGLA[!is.finite(simulationResults$BGoverGLA)] <- NA

simulationResults <- simulationResults[simulationResults$cloneSubtype != "NA_NA_NA",]

# Adds a column indicating to which replicate number the clone belongs
simulationResults$Replicate <- 1
for (i in c(1:replicateNumber)){
  simulationResults$Replicate[(1 + (i-1)*numberClonesPerReplicate) : (numberClonesPerReplicate + (i-1) * numberClonesPerReplicate)] <- i
}

# 9. Compare the simulated clones to the observed ones --------------------------------------------------------------

pdf(file = paste(simulationName,"_ComparisonWithObservedLineages.pdf", sep = ""), width = 8, height = 11)
par(mar=c(5,4,2,2)+0.1, mfrow = c(4,3))

# Plots the input parameters
plot(generationpBG ~ generations, xlim = c(0,25),  xlab = myXlab, main = "Multipoltent progenitor: \nProbabilities of differentiating into the different types", cex.main = 0.8,
     type = 'o', ylim = c(0,1), col = myColors[1], las = 1, ylab = "Probability of producing each astrocyte type", pch = 16)
points(generationpGLA ~ generations,  type = 'o', col = myColors[2], pch = 16)
points(generationpWMA ~ generations,  type = 'o', col = myColors[3], pch = 16)


points(y = firstPostMitoticCell$probability, x= rep(0, length(firstPostMitoticCell$probability)), col = myColors, pch =17)

abline(v = P0Gen , col = "black", lty = 2, lwd = 2)
legend("right", legend = c("WMA", "BG", "Gla"), lty = 1, col = c(myColors[3], myColors[1], myColors[2]),
       bty = "n", cex = 0.8, lwd =2)
legend(x = 14, y = 0.42, legend = "Probabilities for \nsingle cell clones", pch = 17, bty = 'n', cex = 0.8)

##### Plot proba prolif
generationpAstro <-  1 - generationpamplifMP

plot(generationpAstro ~ generations, xlim = c(0,25),  xlab = myXlab, cex.main = 0.9, main = "Cell cycle exit",
     type = 'o', ylim = c(0,1), col = "black", las = 1, ylab = "Probability of cell cycle exit", pch = 16)

if(progType == "MP_WMAP"){
  points(generationpAstro.WMAP ~ generations, type = 'o', col = "black", pch = 18, lty = 2)
  text(x = 5, y = 0.03, labels = paste("Freqency MP: ", freqMP,",\nFrequency WMAP: ", freqWMA.P ), cex = 0.8)
  legend("bottomright", legend = c("MP","WMAP"), lty = c(1,2), bty = 'n')
}

plot.new()

####
## Clone Size
####

myBreaks <- c(0,5,10,15,20, 50, 100, 150, max(observedData$clone_size, simulationResults$CloneSize)+1)
#myBreaks <- seq(1,max(observedData$clone_size, simulationResults$CloneSize)+1,2)
myBreaks <- sort(myBreaks)

observedData$clone_size_Bin <- .bincode(observedData$clone_size, breaks = myBreaks, right = TRUE, include.lowest = T)
simulationResults$Clone_size_Bin <- .bincode(simulationResults$CloneSize, breaks = myBreaks, right = TRUE, include.lowest = T)

observedData$clone_size_Bin <- factor(observedData$clone_size_Bin,levels = seq(1,length(myBreaks)-1,1))
simulationResults$Clone_size_Bin <- factor(simulationResults$Clone_size_Bin,levels = seq(1,length(myBreaks)-1,1))

replicateCloneSizeBin <- meanSDReplicate(simulationResults, variable = "Clone_size_Bin")

propTab <- table(observedData$clone_size_Bin)/length(observedData$clone_size_Bin)
propTab <- rbind(replicateCloneSizeBin[[1]], propTab)
names(propTab)  <- myBreaks[-1]

mp <- barplot(propTab, xlab = "Binned clone size", las = 2, main = "Clone size at P30\n highest value of each bin is shown", cex.main = 0.9,
              cex.axis = 1, cex.names = 0.8, ylab = "proportions of clones", beside = T, names.arg =  myBreaks[-1], col = c("gray", "white"), ylim = c(0, max(propTab) + max(replicateCloneSizeBin[[2]])))
error.bar(x = mp[1,], y = replicateCloneSizeBin[[1]], upper = replicateCloneSizeBin[[2]],length = 0.02)
legend("topright", legend = c("simulated", "observed"), fill = c("gray", "white"), bty = "n")


myError <- sum((propTab[1,]-propTab[2,])^2)

####
### Proportions clone types
####

replicateCloneType <- meanSDReplicate(simulationResults, variable = "cloneType")

mp <- barplot(rbind(replicateCloneType[[1]], prop.table(table(observedData$clone_type_NoCNA))), beside = T, ylim = c(0, 0.7),
              ylab = "Proportion of clones", col = c("gray", "white"), main = "Clone type P30", las = 1, cex.main = 0.9)
error.bar(mp[1,], y = replicateCloneType[[1]], upper = replicateCloneType[[2]], length = 0.02)

# Stats:
tab <- rbind(table(simulationResults$cloneType), table(observedData$clone_type))
chisq.test(t(tab))


####
### Proportions clone subtypes
####


replicateCloneSubtype <- meanSDReplicate(simulationResults, variable = "cloneSubtype")
propTab <- rbind(replicateCloneSubtype[[1]], prop.table(table(observedData$clone_subtype_NoCNA)))
mp <- barplot(rbind(replicateCloneSubtype[[1]], prop.table(table(observedData$clone_subtype_NoCNA))),  ylim = c(0,0.6), cex.main = 0.9,
              beside = T, ylab = "proportions of clones", col = c("gray", "white"), main = "Clone subtype P30", cex.names = 0.7, las = 2)
error.bar(mp[1,], y = replicateCloneSubtype[[1]], upper = replicateCloneSubtype[[2]], length = 0.02)

# Stats:
tab <- rbind(table(simulationResults$cloneSubtype), table(observedData$clone_subtype_NoCNA))
chisq.test(t(tab))

myError <- sum((propTab[1,]-propTab[2,])^2)

####
## Proportions astrocytes
####

## Overall
totalAstroSim <- with(simulationResults, tapply(simulationResults$CloneSize, Replicate, FUN = sum, na.rm = T))

pWMASim <- with(simulationResults, tapply(WMA, INDEX = Replicate, FUN = sum, na.rm = T)) /totalAstroSim
pBGSim <- with(simulationResults, tapply(BG, INDEX = Replicate, FUN = sum, na.rm = T)) /totalAstroSim
pGLASim <- with(simulationResults, tapply(GLA, INDEX = Replicate, FUN = sum, na.rm = T)) /totalAstroSim

meanpWMASim <- mean(pWMASim)
meanpBGSim  <- mean(pBGSim)
meanpGLASim <- mean(pGLASim) 
sdpWMASim <- sd(pWMASim)
sdpBGSim  <- sd(pBGSim)
sdpGLASim <- sd(pGLASim) 

pWMAObs <- with(observedData, sum(number_of_WMA, na.rm = T)/totalAstro)
pBGObs <- with(observedData, sum(number_of_BG, na.rm = T)/totalAstro)
pGLAObs <- with(observedData, sum(number_of_GLA, na.rm = T)/totalAstro)

propTab <- matrix(data = c(meanpBGSim, meanpGLASim, meanpWMASim, pBGObs, pGLAObs, pWMAObs), nrow = 3, ncol = 2)
mp <- barplot(matrix(data = c(meanpBGSim, meanpGLASim, meanpWMASim, pBGObs, pGLAObs, pWMAObs), nrow = 3, ncol = 2), col = myColors,
              names = c("simulated","observed"), ylab = "proportion of astrocytes", las = 1,
              main = "Overall proportion of \n astrocyte types P30", cex.main = 0.9)
error.bar(x = rep(mp[1], times = 3), y = cumsum(c(meanpBGSim, meanpGLASim, meanpWMASim)), upper = c(sdpWMASim,sdpBGSim, sdpGLASim))

myError <- sum((propTab[,1]-propTab[,2])^2)


## Analysis per clone
simTripleClones <- simulationResults$cloneSubtype == "BG_GLA_WMA"
obsTripleClones <- observedData$clone_subtype == "BG_GLA_WMA"

ratioAstroMean <- matrix(data = c(exp(mean(log(simulationResults$GLAoverWMA[simTripleClones]))),
                                  exp(mean(log(simulationResults$BGoverWMA[simTripleClones]))),
                                  exp(mean(log(simulationResults$BGoverGLA[simTripleClones]))),
                                  exp(mean(log(observedData$GLAoverWMA[obsTripleClones]))),
                                  exp(mean(log(observedData$BGoverWMA[obsTripleClones]))),
                                  exp(mean(log(observedData$BGoverGLA[obsTripleClones])))), nrow = 3, ncol = 2)

ratioAstroSEM <- matrix(data = c(exp(sd(log(simulationResults$GLAoverWMA[simTripleClones])/sqrt(sum(simTripleClones)))),
                                 exp(sd(log(simulationResults$BGoverWMA[simTripleClones])/sqrt(sum(simTripleClones)))),
                                 exp(sd(log(simulationResults$BGoverGLA[simTripleClones])/sqrt(sum(simTripleClones)))),
                                 exp(sd(log(observedData$GLAoverWMA[obsTripleClones])/sqrt(sum(obsTripleClones)))),
                                 exp(sd(log(observedData$BGoverWMA[obsTripleClones])/sqrt(sum(obsTripleClones)))),
                                 exp(sd(log(observedData$BGoverGLA[obsTripleClones])/sqrt(sum(obsTripleClones))))), nrow = 3, ncol = 2)

ratioAstroSD <- matrix(data = c(exp(sd(log(simulationResults$GLAoverWMA[simTripleClones]))),
                                exp(sd(log(simulationResults$BGoverWMA[simTripleClones]))),
                                exp(sd(log(simulationResults$BGoverGLA[simTripleClones]))),
                                exp(sd(log(observedData$GLAoverWMA[obsTripleClones]))),
                                exp(sd(log(observedData$BGoverWMA[obsTripleClones]))),
                                exp(sd(log(observedData$BGoverGLA[obsTripleClones])))), nrow = 3, ncol = 2)


#boxplot(simulationResults$WMAoverGLA,observedData$WMAoverGLA, simulationResults$WMAoverBG, observedData$WMAoverBG, simulationResults$GLAoverBG,observedData$GLAoverBG)
myPlot <- barplot(t(ratioAstroMean), beside = T, names = c("GLAoverWMA","BGoverWMA", "BGoveGLA"), col = c("gray", "white"), ylab = "Proportion of astrocyte types per clone (geometric mean)", 
                  main = "Astrocyte type ratios \n triple clones P30", las = 2, cex.names = 0.7, cex.lab = 0.8, ylim = c(0, 8), cex.main = 0.9)

error.bar(x = myPlot, y = t(ratioAstroMean), upper = t(ratioAstroSD))


####
# Comparison P0
####
if(age == "E12"){
  observedDataP0 <- Clones_E12_P0
}
if(age == "E14"){
  observedDataP0 <- Clones_E14_P0
}

table(observedDataP0$clone_subtype.1)
table(observedDataP0$clone_type)


# Proportions of heterogeneous clones in the simulations at generation 6:
storeCloneTableG6 <- storeCloneTable[storeCloneTable$timepoint <= 6 & storeCloneTable$type != "MP", ] # takes only prenatal time points and the differentiated astrocytes
storeCloneTableG6$type <- factor(storeCloneTableG6$type, levels = c("BG", "GLA","WMA"))
simulationResultsG6 <- t(matrix(with(storeCloneTableG6, unlist(tapply(type, INDEX = cloneID, FUN = table))), nrow = 3, ncol =length(unique(storeCloneTableG6$cloneID))))
simulationResultsG6 <- as.data.frame(simulationResultsG6)
names(simulationResultsG6) <- c("BG", "GLA","WMA")

simulationResultsG6$CloneSize <- rowSums(simulationResultsG6)

simulationResultsG6$cloneSubtype <- apply(simulationResultsG6[,1:3], 1, FUN = getCloneSubtype)
simulationResultsG6$cloneType <- apply(simulationResultsG6[,1:3], 1, FUN = getCloneType)


simulationResultsG6$GLAoverWMA <- simulationResultsG6$GLA / simulationResultsG6$WMA
simulationResultsG6$GLAoverWMA[!is.finite(simulationResultsG6$GLAoverWMA)] <- NA

simulationResultsG6$BGoverWMA <- simulationResultsG6$BG / simulationResultsG6$WMA
simulationResultsG6$BGoverWMA[!is.finite(simulationResultsG6$BGoverWMA)] <- NA

simulationResultsG6$BGoverGLA <- simulationResultsG6$BG / simulationResultsG6$GLA
simulationResultsG6$BGoverGLA[!is.finite(simulationResultsG6$BGoverGLA)] <- NA

simulationResultsG6$Replicate <- 1
for ( i in c(1:replicateNumber)){
  simulationResultsG6$Replicate[(1 + (i-1)*numberClonesPerReplicate) : (numberClonesPerReplicate + (i-1) * numberClonesPerReplicate)] <- i
}

replicateCloneTypesP0 <- meanSDReplicate(simulationResultsG6, variable = "cloneType")

plot.new()
mp <- barplot(rbind(replicateCloneTypesP0[[1]], prop.table(table(observedDataP0$clone_type))), 
              beside = T, ylab = "Proportion of clones", col = c("gray", "white"), main = "Clone type P0", las = 1, cex.main = 0.9)
error.bar(mp[1,],y = replicateCloneTypesP0[[1]], upper = replicateCloneTypesP0[[2]])
#legend("topright", legend = c("simulated", "observed"), fill = c("gray", "white"), bty = "n")

tab <- rbind(table(simulationResultsG6$cloneType), table(observedDataP0$clone_type))
chisq.test(t(tab))



# Proportions of cortical vs PWM clones at generation 6. (VZ clones ignored)
summary(as.factor(simulationResultsG6$cloneSubtype))

simulationResultsG6$cloneSubtype2 <- simulationResultsG6$cloneSubtype
for(i in seq_along(simulationResultsG6$cloneSubtype)){
  if(simulationResultsG6$cloneSubtype[i] %in% c("BG", "BG_GLA", "GLA")){simulationResultsG6$cloneSubtype2[i] <- "Cortical"}
  if(simulationResultsG6$cloneSubtype[i] == "WMA"){simulationResultsG6$cloneSubtype2[i] <- "PWM"}
  if(simulationResultsG6$cloneSubtype[i] %in% c("BG_GLA_WMA", "GLA_WMA", "BG_WMA")){simulationResultsG6$cloneSubtype2[i] <- "Cortical_PWM"}
  
}

replicateCloneSubTypesP0 <- meanSDReplicate(simulationResultsG6, variable = "cloneSubtype2")

observedP0cloneSubtypes <- observedDataP0$clone_subtype[observedDataP0$clone_subtype != "VZ"]
observedP0cloneSubtypes <- droplevels(observedP0cloneSubtypes)
mp <- barplot(rbind(replicateCloneSubTypesP0[[1]], prop.table(table(observedP0cloneSubtypes))), ylim = c(0, 0.6), cex.main = 0.9,
              beside = T, ylab = "Proportions of clones", col = c("gray", "white"), main = "Clone subtypes P0", las = 2, cex.names = 0.7)
error.bar(mp[1,], y = replicateCloneSubTypesP0[[1]], upper = replicateCloneSubTypesP0[[2]])
#legend("topleft", legend = c("simulated", "observed"), fill = c("gray", "white"), bty = "n")


tab <- rbind(table(simulationResultsG6$cloneSubtype2), table(observedP0cloneSubtypes))
chisq.test(t(tab))

dev.off()

# Saves all the data tables
save.image(paste(simulationName,".RData", sep = ""))
