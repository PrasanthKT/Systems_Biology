library(WGCNA)
library(tidyverse)

exprData <- read.csv("HW03_expression.csv", check.names = FALSE)

exprRowNames <- make.unique(exprData[, 1])

exprData <- exprData[, -1]
rownames(exprData) <- exprRowNames

anyDuplicated(rownames(exprData)) 

datExpr <- t(exprData)

traitData <- read.csv("HW03_traits.csv", row.names = 1, check.names = FALSE)

commonSamples <- intersect(rownames(datExpr), rownames(traitData))

datExpr <- datExpr[commonSamples, ]
datTraits <- traitData[commonSamples, ]

dim(datExpr)
dim(datTraits)
head(datTraits)

disableWGCNAThreads()

## Question-1

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot the results
par(mfrow = c(1, 2))

plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit (signed R²)", 
     type = "n",
     main = "Scale Independence")
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], 
     labels = powers, 
     cex = 0.8, col = "red")
abline(h = 0.90, col = "red")


plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels = powers, 
     cex = 0.8, col = "red")



# Question - 2

softPower <- 12


adjacency <- adjacency(datExpr, power = softPower)


TOM <- TOMsimilarity(adjacency)

dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot the gene dendrogram
par(mfrow = c(1,1))
par(cex = 0.6)  
par(mar = c(5, 4, 3, 1)) 

plot(geneTree,
     xlab = "", sub = "",
     main = "Gene Clustering Dendrogram Using TOM",
     labels = FALSE,
     hang = 0.03)



# Question - 3

minModuleSize <- 30 

dynamicModules <- cutreeDynamic(dendro = geneTree,
                                distM = dissTOM,
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicModules)
table(dynamicColors) 


plotDendroAndColors(geneTree, 
                    dynamicColors,
                    groupLabels = "Initial Modules",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Initial Module Colors")


MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

plot(METree, 
     main = "Clustering of Module Eigengenes", 
     xlab = "", 
     sub = "")
abline(h = 0.25, col = "red")  # Merge threshold (you can adjust)


merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)


mergedColors <- merge$colors
mergedMEs <- merge$newMEs


plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors),
                    groupLabels = c("Original Modules", "Merged Modules"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Dendrogram with Merged Module Colors")


moduleColors <- mergedColors
MEs <- mergedMEs


# Question - 4

nSamples <- nrow(datExpr)

MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCor
moduleTraitPvalue

# Question - 5

textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)

dim(textMatrix) <- dim(moduleTraitCor)

# Plot the heatmap
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = colnames(datTraits),
  yLabels       = colnames(MEs),
  ySymbols      = colnames(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.9,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships"
)

