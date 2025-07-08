# Organ Volume Analysis
# Rachel Roston, Ph.D.

library(ANTsR)
library(ggplot2)
library(ggsignif)
library(dplyr)

# Load ScanID & Genotype table
metadata = read.csv("./Data/scan_genotypes.csv")
metadata$Genotype = relevel(as.factor(metadata$Genotype), "WT/WT")

# Load label lookup table
labelnames = read.table("./Data/Labels-15.ctbl")[,1:2]
labelnames = rbind(c(0, "background"), labelnames)
names(labelnames) <- c("Label", "Organ")

# Load Images and Labels
img.dir = "./Data/Volumes/"
label.dir = "./Data/Labels/"
mask.dir = "./Data/Masks/"

image.paths = unlist(lapply(metadata$ScanID, function(X){ dir(img.dir, pattern = X, full.names = T)}))
label.paths = unlist(lapply(metadata$ScanID, function(X){ dir(label.dir, pattern = X, full.names = T)}))
mask.paths = unlist(lapply(metadata$ScanID, function(X){ dir(mask.dir, pattern = X, full.names = T)}))

images = lapply(image.paths, antsImageRead)
labels = lapply(label.paths, antsImageRead)
masks = lapply(mask.paths, antsImageRead)

# Make table of organ volumes
labelVolumes = NULL
for(i in 1:length(images)){
  tmp= cbind(metadata[i, c("ScanID", "Genotype")], labelnames[2:16, "Organ"], labelStats(images[[i]], labels[[i]])[2:16,"Volume"])
  labelVolumes = rbind(labelVolumes, tmp)
}
colnames(labelVolumes) <- c("ScanID", "Genotype", "OrganLabel", "Volume")

maskVolumes = data.frame()
for(i in 1:length(images)){
  tmp = cbind.data.frame( ScanID = metadata$ScanID[i], Genotype = metadata$Genotype[i], Volume = labelStats(images[[i]], masks[[i]])[1,"Volume"])
  maskVolumes = rbind(maskVolumes, tmp)
}

labelVolList = split(labelVolumes, labelVolumes$OrganLabel)
labelVolList = lapply(labelVolList, data.frame, Embryo_Volume = maskVolumes$Volume)
allometryTable = bind_rows(labelVolList)
allometryTable$log_Volume = log(allometryTable$Volume)
allometryTable$log_Embryo_Volume = log(allometryTable$Embryo_Volume)

# DO WT & KO SHOW DIFFERENCE IN EMBRYO SIZE?
## Violin plot of embryo volume
embryoVolPlot <- ggplot(maskVolumes, aes(Genotype, Volume)) + 
                    geom_violin() + 
                    theme_bw()
embryoVolPlot

## F-test
var.test(Volume ~ Genotype, data = maskVolumes)

## T-test
t.test(Volume ~ Genotype, data = maskVolumes)


# WHICH ORGANS HAVE SIGNIFICATLY DIFFERENT VOLUMES IN WT VS KO?
## Violin plot of absolute organ volume
organPlot <- ggplot(labelVolumes, aes(x = Genotype, y = Volume)) +
                geom_violin() +
                theme_bw() +
                theme(legend.position = "none") +
                geom_signif(
                  comparisons = list(c("WT/WT", "MUT/MUT")),
                  test = "t.test",
                  map_signif_level = TRUE, 
                  vjust = 1.8,
                  margin_top = 0.05
                ) +
                facet_wrap(~OrganLabel, nrow = 3, ncol = 5, scales = "free", strip.position = "bottom", shrink = F)
organPlot

## t-test: absolute organ volumes
ttest_AbsOrganVol = data.frame()
for(i in 1:length(labelVolList)){
  ttest_AbsOrganVol[i,"OrganLabel"] = names(labelVolList)[i]
  ttest_AbsOrganVol[i, "p.value"] = signif(t.test(Volume ~ Genotype, data = labelVolList[[i]])$p.value, 3)
  ttest_AbsOrganVol[i, "p<0.05"] = ttest_AbsOrganVol[i, "p.value"] < 0.05
}
ttest_AbsOrganVol

## ANCOVA: difference controlling for body size
relVol_ANCOVA_pvals = data.frame(matrix(nrow = length(labelVolList), ncol = 3))
for(i in 1:length(labelVolList)){
  my_model = aov(formula = log(Volume) ~ Genotype + log(Embryo_Volume), 
                 data = labelVolList[[i]])
  my_summary = summary(my_model)
  relVol_ANCOVA_pvals[i, 1] = names(labelVolList)[i]
  relVol_ANCOVA_pvals[i, 2:3] = signif(my_summary[[1]][1:2, "Pr(>F)"], 3)
}
colnames(relVol_ANCOVA_pvals) <- c("Organ", rownames(my_summary[[1]])[1:2])
relVol_ANCOVA_pvals

## Single Regression: which organs correlate w/ body size
relVol_singleReg_pvals = data.frame(matrix(nrow = length(labelVolList), ncol = 2))
for(i in 1:length(labelVolList)){
  my_model = glm(formula = log(Volume) ~ log(Embryo_Volume), 
                 data = labelVolList[[i]])
  my_summary = summary(my_model)
  relVol_singleReg_pvals[i, 1] = names(labelVolList)[i]
  relVol_singleReg_pvals[i, 2] = signif(my_summary$coefficients[2, "Pr(>|t|)"], 3)
}
colnames(relVol_singleReg_pvals) <- c("Organ", "p.val")

## Multiple Regression (likely not relevant to TBM)
relVol_multiReg_pvals = data.frame(matrix(nrow = length(labelVolList), ncol = 4))
for(i in 1:length(labelVolList)){
  my_model = glm(formula = log(Volume) ~ Genotype + log(Embryo_Volume), 
                 data = labelVolList[[i]])
  my_summary = summary(my_model)
  relVol_multiReg_pvals[i, 1] = names(labelVolList)[i]
  relVol_multiReg_pvals[i, 2:4] = signif(my_summary$coefficients[, "Pr(>|t|)"], 3)
}
colnames(relVol_multReg_pvals) <- c("Organ", names(my_summary$coefficients[, "Pr(>|t|)"]))


## Organ & body size scatter plot
evp <- ggplot(data = allometryTable, aes(x = log_Embryo_Volume, y = log_Volume, shape = Genotype, colour = Genotype)) +
  geom_point() +
  scale_color_grey(start = 0.2, end = 0.6) +
  theme_bw() +
  facet_wrap(~OrganLabel, nrow = 3, ncol = 5, scales = "free")
evp
# evp +  
#   geom_smooth(method = lm,se=F,aes(group=1)) +
#   geom_smooth(method = lm,se=F)
