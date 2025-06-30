# Registration intensityOnly
# Rachel Roston

experiment = "labelIntensity"

library(ANTsR)
  out.dir = paste0("./Results/", experiment, "/Transforms/")
  out.dir.subfolders = paste0(out.dir, c("WT", "KO"))
  sapply(X = out.dir.subfolders, FUN = dir.create, recursive = T)

  metadata = read.csv("./Data/scan_genotypes.csv")
  
  img.dir = "./Data/Volumes/"
  label.dir = "./Data/Labels/"
  mask.dir = "./Data/Masks/"
  
  ref.image = antsImageRead("./Data/Template/maskedtemplate0__lowRes.nii.gz")
  ref.label = antsImageRead("./Data/Template/maskedtemplate0__lowRes-15-label.nii.gz")
  ref.mask = antsImageRead("./Data/Template/maskedtemplate0__lowRes-wholebody-label.nii.gz")
  
  for(i in 1:nrow(metadata)){
    ScanID = metadata$ScanID[i]
    out.prefix = paste0(out.dir, "/", metadata$Subdirectory[i], "/", ScanID, "__", experiment, "_")
    
    mov.image = antsImageRead(paste0(img.dir, ScanID, "__rec-18um-tx-low.nrrd"))
    mov.label = antsImageRead(paste0(label.dir, ScanID, "__tx-low-15-label.nii.gz"))
    mov.mask = antsImageRead(paste0(mask.dir, ScanID, "__tx-low-wholebody-label.nii.gz"))
    
    # Step1: intensity only
      # linear = antsRegistration(fixed = ref.image,
      #                           moving = mov.image,
      #                           typeofTransform = "Similarity",
      #                           mask = ref.mask,
      #                           movingMask = mov.mask)
      # step1 = antsRegistration(fixed = ref.image,
      #                           moving = mov.image,
      #                           typeofTransform = "antsRegistrationSyN[so]",
      #                           mask = ref.mask,
      #                           movingMask = mov.mask,
      #                           initialTransforms = linear$fwdtransforms,
      #                           outprefix = paste0(out.prefix, "Step1_"))
      
    # Step1: label-intensity
      step1 = labelImageRegistration(fixedLabelImages = ref.label,
                                     movingLabelImages = mov.label,
                                     fixedIntensityImages = ref.image,
                                     movingIntensityImages = mov.image,
                                     fixedMask = ref.mask,
                                     movingMask = mov.mask,
                                     initialTransforms = "similarity",
                                     typeOfDeformableTransform = "antsRegistrationSyN[so]",
                                     outputPrefix = paste0(out.prefix, "Step1_"))
      
    tx.img.step1Linear = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step1$fwdtransforms[2],
                                       interpolator = "linear")
    antsImageWrite(tx.img.step1Linear, paste0(out.prefix, "Step1_linearTx_fwdWarpedimg.nrrd"))
    
    tx.img.step1 = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step1$fwdtransforms,
                                       interpolator = "linear")
    antsImageWrite(tx.img.step1, paste0(out.prefix, "Step1_Warp1Tx_fwdWarpedimg.nrrd"))
    
    # Step 2
    step2 = antsRegistration(fixed = ref.image,
                              moving = mov.image,
                              typeofTransform = "antsRegistrationSyN[so]",
                              mask = ref.mask,
                              movingMask = mov.mask,
                              initialTransforms = step1$fwdtransforms,
                              outprefix = paste0(out.prefix, "Step2_"))

    tx.img.step2 = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step2$fwdtransforms,
                                       interpolator = "linear")
    antsImageWrite(tx.img.step2, paste0(out.prefix, "Step2_Warp2Tx_fwdWarpedImg.nrrd"))
    
    
    antsApplyTransforms(fixed= ref.image, 
                        moving = mov.image, 
                        transformlist = step2$fwdtransforms, 
                        compose = paste0(out.prefix, "combinedFwdWarps-"))
    
    # Inverse Labels
    inverse.labels = antsApplyTransforms(fixed = mov.image, 
                                         moving = ref.label,
                                         transformlist = c(step1$invtransforms[1],
                                                           step1$invtransforms[2],
                                                           step2$invtransforms[2]),
                                         whichtoinvert = c(T,F,F), 
                                         interpolator = "genericLabel")
    antsImageWrite(inverse.labels, paste0(out.prefix, "inverse-label.nii.gz"))
    
    remove(mov.image, mov.label, mov.mask, 
           step1, step2, 
           tx.img.step1Linear, tx.img.step1, tx.img.step2,
           inverse.labels)
    gc()
  }
  