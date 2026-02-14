

library(ANTsR)
library(ggplot2)
library(ggsignif)
library(dplyr)

img.dir = "./Data/Volumes/"
mask.dir = "./Data/Masks/"

ma.label.dir = "./Data/Labels/"
io.label.dir = "./Results/intensityOnly/invTx_Labels/"
li.label.dir = "./Results/labelIntensity/invTx_Labels/"

# Load label lookup table
labelnames = read.table("./Data/Labels-15.ctbl")[,1:2]
labelnames = rbind(c(0, "background"), labelnames)
names(labelnames) <- c("Label", "Organ")
labelnames$OrganName = gsub("_left", " (left)", labelnames$Organ)
labelnames$OrganName = gsub("_right", " (right)", labelnames$OrganName)

metadata = read.csv("./Data/scan_genotypes.csv")
ko = which(metadata$Genotype == "MUT/MUT")
wt = which(metadata$Genotype == "WT/WT")
metadata = metadata[c(wt, ko), ]

imgs = lapply(metadata$ScanID, function(X){antsImageRead(dir(img.dir, pattern = X, full.names = T))})
masks = lapply(metadata$ScanID, function(X){antsImageRead(dir(mask.dir, pattern = X, full.names = T))})
ma.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(ma.label.dir, pattern = X, full.names = T, recursive = T))})
io.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(io.label.dir, pattern = paste0(X,".*step2"), full.names = T, recursive = T))})
li.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(li.label.dir, pattern = paste0(X,".*step2"), full.names = T, recursive = T))})


ma.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], labelStats(imgs[[X]], ma.labels[[X]]))))
io.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], labelStats(imgs[[X]], io.labels[[X]]))))
li.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], labelStats(imgs[[X]], li.labels[[X]]))))

volumes = data.frame(LabelType = rep(c("m", "i", "l"), each = nrow(ma.vols)),
                     rbind(ma.vols, io.vols, li.vols))
volumes$GenotypeLabel = paste0(volumes$Subdirectory, "-", volumes$LabelType)




organVolPlot <- ggplot(volumes, aes(x = GenotypeLabel, y = Volume)) +
             geom_boxplot() +
             theme_bw() +
             theme(legend.position = "none") +
             facet_wrap(~LabelValue, nrow = 4, ncol = 4, scales = "free", strip.position = "bottom", shrink = F) +
             theme(text = element_text(size = 4))
ggsave("./Results/organVols_violinPlot_2.tif", organVolPlot)


io.overlap = do.call(rbind, lapply(1:length(ma.labels), function(X) data.frame(metadata[X,], RegType = "Intensity Only", labelOverlapMeasures(ma.labels[[X]], io.labels[[X]]))))
li.overlap = do.call(rbind, lapply(1:length(ma.labels), function(X) data.frame(metadata[X,], RegType = "Label-Informed", labelOverlapMeasures(ma.labels[[X]], li.labels[[X]]))))

diceCoeff.table = rbind(io.overlap, li.overlap)
diceCoeff.table$Subdirectory <- relevel(factor(diceCoeff.table$Subdirectory), ref = "WT")
diceCoeff.table$GenotypeRegtype <- factor(paste0(diceCoeff.table$Subdirectory, "\n", diceCoeff.table$RegType),
                                          levels = c("WT\nIO", "KO\nIO", "WT\nLI", "KO\nLI"), 
                                          ordered = T)
for(i in 1:nrow(diceCoeff.table)){
  if(diceCoeff.table$Label[i] == "All") {
    diceCoeff.table$OrganName[i] <- "All"
  } else {
    diceCoeff.table$OrganName[i] <- labelnames[which(labelnames$Label == diceCoeff.table$Label[i]), "OrganName"]
  }
}
diceCoeff.table$OrganName <- relevel(factor(diceCoeff.table$OrganName), ref = "All")


DicePlot1<-ggplot(diceCoeff.table, aes(x = GenotypeRegtype, y = MeanOverlap)) +
  geom_boxplot() +
  ylim(0,1) +
  labs(x = "", y = "Dice") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~OrganName, nrow = 4, ncol = 4, scales = "fixed", 
             strip.position = "top", 
             shrink = F, 
             axes = "margins", 
             axis.labels = "margins") +
  theme(text = element_text(size = 9)) +
  theme(axis.text.x = element_text(hjust = 0.5))
ggsave("./Results/dicePlot_test.tif", DicePlot1)


DicePlot <- ggplot(diceCoeff.table, aes(x = OrganName, y = MeanOverlap, color = OrganName)) +
  geom_boxplot() +
  scale_color_grey() +
  theme_bw() +
  labs(x = "", y = "Dice") +
  facet_grid(cols = vars(Subdirectory), rows = vars(RegType), switch = "both") +
  theme(legend.position = "right",
        text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        axis.text.x = element_blank(),      # Remove x-axis text
        axis.ticks.x = element_blank()) 
ggsave("./Results/dice_test2.tif", DicePlot)
