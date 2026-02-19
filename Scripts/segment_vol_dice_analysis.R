

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
labelnames$OrganName = gsub("_", " ", labelnames$OrganName)

metadata = read.csv("./Data/scan_genotypes.csv")
ko = which(metadata$Genotype == "MUT/MUT")
wt = which(metadata$Genotype == "WT/WT")
metadata = metadata[c(wt, ko), ]

imgs = lapply(metadata$ScanID, function(X){antsImageRead(dir(img.dir, pattern = X, full.names = T))})
masks = lapply(metadata$ScanID, function(X){antsImageRead(dir(mask.dir, pattern = X, full.names = T))})
ma.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(ma.label.dir, pattern = X, full.names = T, recursive = T))})
io.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(io.label.dir, pattern = paste0(X,".*step2"), full.names = T, recursive = T))})
li.labels = lapply(metadata$ScanID, function(X){antsImageRead(dir(li.label.dir, pattern = paste0(X,".*step2"), full.names = T, recursive = T))})

# Dice


io.overlap = do.call(rbind, lapply(1:length(ma.labels), function(X) data.frame(metadata[X,], RegType = "Intensity Only", labelOverlapMeasures(ma.labels[[X]], io.labels[[X]]))))
li.overlap = do.call(rbind, lapply(1:length(ma.labels), function(X) data.frame(metadata[X,], RegType = "Label-Informed", labelOverlapMeasures(ma.labels[[X]], li.labels[[X]]))))

diceCoeff.table = rbind(io.overlap, li.overlap)
diceCoeff.table$Subdirectory <- relevel(factor(diceCoeff.table$Subdirectory), ref = "WT")
for(i in 1:nrow(diceCoeff.table)){
  if(diceCoeff.table$Label[i] == "All") {
    diceCoeff.table$OrganName[i] <- "All"
  } else {
    diceCoeff.table$OrganName[i] <- labelnames[which(labelnames$Label == diceCoeff.table$Label[i]), "OrganName"]
  }
}
diceCoeff.table$OrganName <- relevel(factor(diceCoeff.table$OrganName), ref = "All")
write.csv(diceCoeff.table, "./Results/dice_table.csv")

DicePlot <- ggplot(diceCoeff.table, aes(x = OrganName, y = MeanOverlap, color = OrganName)) +
  theme_bw() +
  geom_boxplot(linewidth = 0.3, color = "black") +
  labs(x = "", y = "Dice") +
  facet_grid(cols = vars(Subdirectory), rows = vars(RegType), switch = "both") +
  theme(legend.position = "none",
        text = element_text(size = 9),
        axis.text.x = element_text(margin=margin(t=0.5),angle = 45, vjust = 1, hjust = 1),      # Remove x-axis text
        #axis.ticks.x = element_blank()
  )
ggsave("./Results/dice_plot.tif", DicePlot, width = 5, height = 4, dpi = 300)


# Volume

ma.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], RegType = "M", labelStats(imgs[[X]], ma.labels[[X]]))))
io.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], RegType = "IO", labelStats(imgs[[X]], io.labels[[X]]))))
li.vols = do.call(rbind, lapply(1:length(imgs), function(X) data.frame(metadata[X,], RegType = "LI", labelStats(imgs[[X]], li.labels[[X]]))))

volumes = rbind(ma.vols, io.vols, li.vols)
for(i in 1:nrow(volumes)){
    volumes$OrganName[i] <- labelnames[which(labelnames$Label == volumes$Label[i]), "OrganName"]
}
volumes = volumes[which(volumes$LabelValue !=0),]
volumes$GenotypeRegtype = factor(paste0(volumes$Subdirectory, "-", volumes$RegType),
                                 levels = c("WT-M", "WT-IO", "WT-LI", "KO-M", "KO-IO", "KO-LI"),
                                 ordered = T)
write.csv(volumes, "./Results/volumes_table.csv")

organVolPlot <- ggplot(volumes, aes(x = GenotypeRegtype, y = Volume)) +
  theme_bw() +
  geom_boxplot(linewidth = 0.3, color = "black") +
  labs(x = "", y = expression("Organ Volume (mm"^3*")")) +
  facet_wrap(~OrganName, nrow = 3, ncol = 5, scales = "free", strip.position = "bottom", shrink = F) +
  theme(legend.position = "none",
        text = element_text(size = 9),
        axis.text.x = element_text(margin=margin(t=0.5),angle = 45, vjust = 1, hjust = 1),      # Remove x-axis text
        #axis.ticks.x = element_blank()
  )
ggsave("./Results/organVols_plot.tif", organVolPlot, width = 7, height = 5, dpi =300)
