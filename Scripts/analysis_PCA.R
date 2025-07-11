# Principal Component Analysis (PCA)
# Rachel Roston, Ph.D.

library(ANTsR)
library(stringr)
library(viridis)

# All forward transforms  
  # intensity-only PCA
  
  experiment = "intensityOnly"
  txPattern = "combinedFwdWarps"
  baseline = "WT"
  
  {
    tx.dir = paste0("./Results/", experiment, "/Transforms")
    
    tx.paths = dir(tx.dir, pattern = txPattern, recursive = T)
    subjects = str_extract(tx.paths, "Scan_....")
    genotypes = str_split_i(tx.paths, "/", 1)
    genotypes = relevel(as.factor(genotypes), ref = baseline)
    
    ref.img = antsImageRead("./Data/Template/maskedtemplate0__lowRes.nii.gz")
    ref.mask = antsImageRead("./Data/Template/maskedtemplate0__lowRes-wholebody-label.nii.gz")
    
    # Output directory
    results.dir = paste0("./Results/", experiment, "/PCA-", txPattern)
    if(! dir.exists(results.dir)) dir.create(results.dir)
    
    # PCA
    pca_mask = iMath(ref.mask, "MD", 15) # dilate to allow PC img deformations
    
    transforms = lapply(paste0(tx.dir, "/", tx.paths), antsImageRead)
    
    pca <- multichannelPCA(x = transforms, 
                           mask = pca_mask, 
                           pcaOption = "randPCA")
    
    # save csv for each PCA output (d, u, v)
    pca_d = data.frame(d = pca$pca$d)
    write.csv(pca_d, paste0(results.dir, "/pca_d.csv"))
    
    pca_u = data.frame(ScanID = subjects, pca$pca$u)
    colnames(pca_u) <- c("ScanID", paste0("PC", 1:(ncol(pca_u)-1)))
    write.csv(pca_u, paste0(results.dir, "/pca_u.csv"))
    
    # save deformed atlas
    PC_warpedImg.dir = paste0(results.dir, "/warpedTemplate")
    if(!dir.exists(PC_warpedImg.dir))dir.create(PC_warpedImg.dir)
    
    for(i in 1:ncol(pca$pca$u)){
      pos.scaled = vectorToMultichannel(500*pca$pca$v[,i], pca_mask)   #this should be the mask you used for your PCA
      pos.scaledTX = antsrTransformFromDisplacementField(pos.scaled)
      pos.warped.img = applyAntsrTransform(pos.scaledTX, 
                                           data=ref.img,
                                           reference= ref.img, 
                                           interpolation = "linear")
      antsImageWrite(pos.warped.img, paste0(PC_warpedImg.dir, "/PC", str_pad(i, width = 2, side = "left", 0), "_", experiment, "_pos-500_warpedTemplate.nrrd"))
      
      neg.scaled = vectorToMultichannel(-500*pca$pca$v[,i], pca_mask)   #this should be the mask you used for your PCA
      neg.scaledTX = antsrTransformFromDisplacementField(neg.scaled)
      neg.warped.img = applyAntsrTransform(neg.scaledTX, 
                                           data=ref.img,
                                           reference= ref.img, 
                                           interpolation = "linear")
      antsImageWrite(neg.warped.img, paste0(PC_warpedImg.dir, "/PC", str_pad(i, width = 2, side = "left", 0), "_", experiment, "_neg-500_warpedTemplate.nrrd"))
      
    }
    
    # save plots
    wt = which(genotypes == "WT")
    ko = which(genotypes == "KO")
    
    ## Scree plot
    par(family="")
    
    png(filename = paste0(results.dir, "/scree.png"),
        #res = 300,
        units = "px",
        height = 480,
        width = 880)
    
    eigen.sum=sum(pca$pca$d)
    PC.percent = pca$pca$d/eigen.sum*100
    plot(x = 1:length(pca$pca$d),
         y = PC.percent,
         xlab = "PC",
         ylab = "% Variance Explained",
         main = paste0(experiment, ", ", txPattern))
    dev.off()
    
    ## PC-PC plot
    maxPC = floor(length(pca$pca$d)/2)
    
    for(p in seq(1, maxPC, by = 2)){
      PCs = c(p, p+1)
      
      par(family="")
      
      png(filename = paste0(results.dir, "/PC_plot", PCs[1], "-", PCs[2], ".png"),
          #res = 300,
          units = "px",
          height = 480,
          width = 480)
      
      par(mar= c(5,5,5,5), xpd = TRUE)
      
      plot(x = pca$pca$u[,PCs[1]], 
           y = pca$pca$u[,PCs[2]], 
           pch=32, 
           xlab = paste0("PC", PCs[1], " (", signif(PC.percent, digits = 3)[PCs[1]], "%)"), 
           ylab = paste0("PC", PCs[2], " (", signif(PC.percent, digits = 3)[PCs[2]], "%)"),
           main = paste0(experiment, ", ", txPattern)
      )
      # plot wildtype
      points(x = pca$pca$u[wt,PCs[1]], 
             y = pca$pca$u[wt,PCs[2]], pch=17, col="black")
      # plot KOs
      points(x = pca$pca$u[ko,PCs[1]], 
             y = pca$pca$u[ko,PCs[2]], pch=16, col="red")
      legend('topright',
             legend = levels(genotypes),
             col = c("black", "red"),
             pch = c(17, 16))
      dev.off()
    }
  }
  
  # label-intensity PCA
  
  experiment = "labelIntensity"
  txPattern = "combinedFwdWarps"
  baseline = "WT"
  
  {
    tx.dir = paste0("./Results/", experiment, "/Transforms")
    
    tx.paths = dir(tx.dir, pattern = txPattern, recursive = T)
    subjects = str_extract(tx.paths, "Scan_....")
    genotypes = str_split_i(tx.paths, "/", 1)
    genotypes = relevel(as.factor(genotypes), ref = baseline)
    
    ref.img = antsImageRead("./Data/Template/maskedtemplate0__lowRes.nii.gz")
    ref.mask = antsImageRead("./Data/Template/maskedtemplate0__lowRes-wholebody-label.nii.gz")
    
    # Output directory
    results.dir = paste0("./Results/", experiment, "/PCA-", txPattern)
    if(! dir.exists(results.dir)) dir.create(results.dir)
    
    # PCA
    pca_mask = iMath(ref.mask, "MD", 15) # dilate to allow PC img deformations
    
    transforms = lapply(paste0(tx.dir, "/", tx.paths), antsImageRead)
    
    pca <- multichannelPCA(x = transforms, 
                           mask = pca_mask, 
                           pcaOption = "randPCA")
    
    # save csv for each PCA output (d, u, v)
    pca_d = data.frame(d = pca$pca$d)
    write.csv(pca_d, paste0(results.dir, "/pca_d.csv"))
    
    pca_u = data.frame(ScanID = subjects, pca$pca$u)
    colnames(pca_u) <- c("ScanID", paste0("PC", 1:(ncol(pca_u)-1)))
    write.csv(pca_u, paste0(results.dir, "/pca_u.csv"))
    
    # save deformed atlas
    PC_warpedImg.dir = paste0(results.dir, "/warpedTemplate")
    if(!dir.exists(PC_warpedImg.dir))dir.create(PC_warpedImg.dir)
    
    for(i in 1:ncol(pca$pca$u)){
      pos.scaled = vectorToMultichannel(500*pca$pca$v[,i], pca_mask)   #this should be the mask you used for your PCA
      pos.scaledTX = antsrTransformFromDisplacementField(pos.scaled)
      pos.warped.img = applyAntsrTransform(pos.scaledTX, 
                                           data=ref.img,
                                           reference= ref.img, 
                                           interpolation = "linear")
      antsImageWrite(pos.warped.img, paste0(PC_warpedImg.dir, "/PC", str_pad(i, width = 2, side = "left", 0), "_", experiment, "_pos-500_warpedTemplate.nrrd"))
      
      neg.scaled = vectorToMultichannel(-500*pca$pca$v[,i], pca_mask)   #this should be the mask you used for your PCA
      neg.scaledTX = antsrTransformFromDisplacementField(neg.scaled)
      neg.warped.img = applyAntsrTransform(neg.scaledTX, 
                                           data=ref.img,
                                           reference= ref.img, 
                                           interpolation = "linear")
      antsImageWrite(neg.warped.img, paste0(PC_warpedImg.dir, "/PC", str_pad(i, width = 2, side = "left", 0), "_", experiment, "_neg-500_warpedTemplate.nrrd"))
      
    }
    
    # save plots
    wt = which(genotypes == "WT")
    ko = which(genotypes == "KO")
    
    ## Scree plot
    par(family="")
    
    png(filename = paste0(results.dir, "/scree.png"),
        #res = 300,
        units = "px",
        height = 480,
        width = 880)
    
    eigen.sum=sum(pca$pca$d)
    PC.percent = pca$pca$d/eigen.sum*100
    plot(x = 1:length(pca$pca$d),
         y = PC.percent,
         xlab = "PC",
         ylab = "% Variance Explained",
         main = paste0(experiment, ", ", txPattern))
    dev.off()
    
    ## PC-PC plot
    maxPC = floor(length(pca$pca$d)/2)
    
    for(p in seq(1, maxPC, by = 2)){
      PCs = c(p, p+1)
      
      par(family="")
      
      png(filename = paste0(results.dir, "/PC_plot", PCs[1], "-", PCs[2], ".png"),
          #res = 300,
          units = "px",
          height = 480,
          width = 480)
      
      par(mar= c(5,5,5,5), xpd = TRUE)
      
      plot(x = pca$pca$u[,PCs[1]], 
           y = pca$pca$u[,PCs[2]], 
           pch=32, 
           xlab = paste0("PC", PCs[1], " (", signif(PC.percent, digits = 3)[PCs[1]], "%)"), 
           ylab = paste0("PC", PCs[2], " (", signif(PC.percent, digits = 3)[PCs[2]], "%)"),
           main = paste0(experiment, ", ", txPattern)
      )
      # plot wildtype
      points(x = pca$pca$u[wt,PCs[1]], 
             y = pca$pca$u[wt,PCs[2]], pch=17, col="black")
      # plot KOs
      points(x = pca$pca$u[ko,PCs[1]], 
             y = pca$pca$u[ko,PCs[2]], pch=16, col="red")
      legend('topright',
             legend = levels(genotypes),
             col = c("black", "red"),
             pch = c(17, 16))
      dev.off()
    }
  }