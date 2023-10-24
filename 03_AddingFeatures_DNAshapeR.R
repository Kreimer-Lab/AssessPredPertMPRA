library(DNAshapeR)
library("Biostrings")
library(dplyr)
library(tibble)

fa <- readDNAStringSet("./locs.fasta")
sequence_names <- names(fa)
pred_inter <- getShape("./locs.fasta", 
                       shapeType = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt' # inter bp
                       ))
pred_intra <- getShape("./locs.fasta", 
                       shapeType = c(
                         'Buckle', 'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', "MGW" #intra bp
                       ))

pred_inter <- 
  lapply(pred_inter, function(X){
    X[is.na(X)] = 0
    rownames(X) = sequence_names
    return(X)
  })
pred_intra <- 
  lapply(pred_intra, function(X){
    X[is.na(X)] = 0
    rownames(X) = sequence_names
    return(X)
  })

pred_inter <- 
  lapply(pred_inter, function(X){
    rowSums(X)/ncol(X)
  })

pred_intra <- 
  lapply(pred_intra, function(X){
    rowSums(X)/ncol(X)
  })


pred_inter <- as.data.frame(pred_inter)
pred_intra <- as.data.frame(pred_intra)

pred <- cbind(pred_inter, pred_intra)

write.csv(pred, "DNAshapeR.csv", quote = F)
