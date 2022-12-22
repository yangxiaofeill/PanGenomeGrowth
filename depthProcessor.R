library(data.table)
library(matrixStats)
library(ggpubr)
getCumuSize <- function(outFile, hapDepth, hapNum, IDIndexes, lineClasses){
  ## get average length
  message("---------- in getCumuSize function -------- ")
  
  message("----- [ get TRUE or FALSE matrix of hapDepth ] -----")
  hapDepthMatrix <- hapDepth > 0
  
  depthClass <- lineClasses
  
  ##################
  
  message("----- [ calculating the cumuSize ] -----")
  hapDepthSize <- data.frame(sample = "", size = 0)
  hapDepthSizeIncrease <- data.frame(sample = "", sizeIncrease = 0, cumuSize = 0)
  hapDepthSize <- hapDepthSize[-1, ]
  hapDepthSizeIncrease <- hapDepthSizeIncrease[-1, ]
  theClass <- c("single", "poly", "common", "core")
  
  currSign <- c()
  currSize <- data.frame(class="", cumuSize = 0)
  currSize <- currSize[-1, ]
  
  for(i in IDIndexes){
    print(i)
    sampleID <- colnames(hapDepth)[i]
    message(paste(sampleID, sep = " "))
    oneHapDep <- hapDepth[, i]
    if(i == IDIndexes[1]){
      currSign <- hapDepthMatrix[, i]
      for(oneClass in theClass){
        message(paste(sampleID, "------", oneClass, sep = " "))
        oneClassSize <- sum(oneHapDep[which(depthClass == oneClass & hapDepthMatrix[, i])])
        hapDepthSizeIncrease <- rbind(hapDepthSizeIncrease,
                                      data.frame(sample = sampleID,
                                                 type = oneClass,
						 sizeIncrease = 0,
                                                 cumuSize = sum(oneClassSize)))
        currSize <- rbind(currSize, data.frame(class = oneClass, cumuSize = sum(oneClassSize)))
      }
    }else{
      newSign <- currSign | hapDepthMatrix[, i]
      oneHapDep <- hapDepth[, i]
      for(oneClass in theClass){
        message(paste(sampleID, "------", oneClass, sep = " "))
        oneClassIncreaseSize <- sum(oneHapDep[which(newSign & !currSign & depthClass == oneClass)])
        currCumuSize <- currSize$cumuSize[which(currSize$class == oneClass)]
        hapDepthSizeIncrease <- rbind(hapDepthSizeIncrease,
                                      data.frame(sample = sampleID,
					         type = oneClass,
                                                 sizeIncrease = oneClassIncreaseSize,
                                                 cumuSize = currCumuSize + oneClassIncreaseSize))
        currSize$cumuSize[which(currSize$class == oneClass)] <- currCumuSize + oneClassIncreaseSize
      }
      currSign <- newSign
    }
  }
  write.table(hapDepthSizeIncrease, file = outFile, row.names = F, col.names=T, quote=F, sep = "\t")
}

args <- commandArgs(trailingOnly = TRUE)
theVCFFile      <- args[1]
sampleNum       <- as.numeric(args[2])
sampleOrderFile <- "xxx"
if(length(args) == 3){
   sampleOrderFile <- args[3]
}

hapNum    <- sampleNum * 2
jarName   <- "PGG"
theV      <- "v1.0"
message(paste("the file is ", theVCFFile, sep = ""))
message("----- begin reading hapDepth file -----")
hapDepth <- fread(paste(theVCFFile, ".hapDepth_", jarName, "_", theV, sep = ""), header = T, stringsAsFactors = F)
lineID <- hapDepth[, "lineIndex"]
lineClass <- hapDepth[, "lineClass"]
hapDepth <- as.matrix(hapDepth[, 1:hapNum])

IDIndexes <- 1:hapNum
if(sampleOrderFile != "xxx"){
   IDindexes <- c()
   theIDs <- fread(sampleOrderFile, header = F, stringsAsFactor = F)
   vcfSampleID <- colnames(hapDepth)
   IDIndexes <- c()
   for(oneID in theIDs$V1){
        oneindex1 <- which(vcfSampleID == paste(oneID, "_H1", sep = ""))
        oneindex2 <- which(vcfSampleID == paste(oneID, "_H2", sep = ""))
        IDIndexes <- c(IDIndexes, oneindex1, oneindex2)
   }
}

message("-------- processing haplotype depth files by depthProcessor.R ------")
message(paste("version is ", theV, sep = ""))
message("Rscript --vanilla depthProcessor.R theVCFFile[required] sampleNum[required] orderedSampleFile[optional]")
outFile <- paste(theVCFFile, ".hapDepth_", jarName, "_", theV, "_cumuSize", sep = "")
getCumuSize(outFile = outFile, hapDepth = hapDepth, hapNum = hapNum, IDIndexes = IDIndexes, lineClasse = lineClass)
message(paste("The output file is ", outFile, sep = ""))
message("The output figure is panGenomeGrowth.pdf, with width = 9, and height = 5")

hapDepthCumu <- fread(outFile, header = T, stringsAsFactors = F)
hapDepthCumu$cumuSize <- round(hapDepthCumu$cumuSize / 1e6, 2)
hapDepthCumu$type <- factor(hapDepthCumu$type, levels = c("single", "poly", "common", "core"))

finalSign <- hapDepthCumu$sample[nrow(hapDepthCumu)]
pdf("panGenomeGrowth.pdf", width = 9, height=5)
ggbarplot(data = hapDepthCumu, x = "sample", y = "cumuSize",
          fill = "type", color = "type",
          title = paste("Pangenome growth, total: ",sum(hapDepthCumu$cumuSize[which(hapDepthCumu$sample == finalSign)]), sep = ""),
          palette = c("#FFC75F","#4D8076", "#C34A36", "#845EC2"),
          lab.col = "black", ylab = "Mbs", xlab = "Haplotypes", x.text = NULL,
          rotate = F)
print(hapDepthCumu[which(hapDepthCumu$sample == finalSign), c("type", "cumuSize")], row.names = F, quote = F)

message("*********THE END ***************")
