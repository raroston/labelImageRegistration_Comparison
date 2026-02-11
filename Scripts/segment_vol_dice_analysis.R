

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
             geom_violin() +
             theme_bw() +
             theme(legend.position = "none") +
             facet_wrap(~LabelValue, nrow = 4, ncol = 4, scales = "free", strip.position = "bottom", shrink = F) +
             theme(text = element_text(size = 4))
ggsave("./Results/organVols_violinPlot.tif", organVolPlot)



