# Registration intensityOnly
# Rachel Roston

library(ANTsR)

# Set Directory for Registration Output
experiment = "intensityOnly"

# Set Up Subfolders for Registration Output
  out.dir = paste0("./Results/", experiment)
  
  out.dir.transforms = paste0(out.dir, "/Transforms")
  out.dir.txVol = paste0(out.dir, "/fwdTx_Volumes")
  out.dir.txLab = paste0(out.dir, "/invTx_Labels")
  sapply(X = paste0(out.dir.transforms, c("/WT", "/KO")), FUN = dir.create, recursive = T)
  sapply(X = paste0(out.dir.txVol, c("/WT", "/KO")), FUN = dir.create, recursive = T)
  sapply(X = paste0(out.dir.txLab, c("/WT", "/KO")), FUN = dir.create, recursive = T)

# Set Input Data Directories for Subjects
  img.dir = "./Data/Volumes/"
  label.dir = "./Data/Labels/"
  mask.dir = "./Data/Masks/"

# Load Input Data
  metadata = read.csv("./Data/scan_genotypes.csv")
  ref.image = antsImageRead("./Data/Template/maskedtemplate0__lowRes.nii.gz")
  ref.label = antsImageRead("./Data/Template/maskedtemplate0__lowRes-15-label.nii.gz")
  ref.mask = antsImageRead("./Data/Template/maskedtemplate0__lowRes-wholebody-label.nii.gz")

# Registration Loop
  for(i in 1:nrow(metadata)){
    ScanID = metadata$ScanID[i]
    out.prefix = paste0(ScanID, "__", experiment, "_")
    
    print(paste0("Starting registration ", i, "/", nrow(metadata), ": ", ScanID))
    start = Sys.time()
    
    mov.image = antsImageRead(paste0(img.dir, ScanID, "__rec-18um-tx-low.nrrd"))
    mov.label = antsImageRead(paste0(label.dir, ScanID, "__tx-low-15-label.nii.gz"))
    mov.mask = antsImageRead(paste0(mask.dir, ScanID, "__tx-low-wholebody-label.nii.gz"))
    
    # Step1
    print(paste("Starting Step1 registration,", Sys.time()))
    linear = antsRegistration(fixed = ref.image,
                              moving = mov.image,
                              typeofTransform = "Similarity",
                              mask = ref.mask,
                              movingMask = mov.mask)
    step1 = antsRegistration(fixed = ref.image,
                              moving = mov.image,
                              typeofTransform = "antsRegistrationSyN[so]",
                              mask = ref.mask,
                              movingMask = mov.mask,
                              initialTransform = linear$fwdtransforms,
                              outprefix = paste0(out.dir.transforms, "/", metadata$Subdirectory[i], "/", out.prefix, "Step1_"))
      
    tx.img.step1Linear = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step1$fwdtransforms[2],
                                       interpolator = "linear")
    antsImageWrite(tx.img.step1Linear, paste0(out.dir.txVol, "/", metadata$Subdirectory[i], "/", out.prefix, "Step1_linearTx_fwdWarpedimg.nrrd"))
    
    tx.img.step1 = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step1$fwdtransforms,
                                       interpolator = "linear")
    antsImageWrite(tx.img.step1, paste0(out.dir.txVol, "/", metadata$Subdirectory[i], "/", out.prefix, "Step1_Warp1Tx_fwdWarpedimg.nrrd"))
    
    # Step 2
    print(paste("Starting Step2 registration,", Sys.time()))
    step2 = antsRegistration(fixed = ref.image,
                              moving = mov.image,
                              typeofTransform = "antsRegistrationSyN[so]",
                              mask = ref.mask,
                              movingMask = mov.mask,
                              initialTransform = step1$fwdtransforms,
                              outprefix = paste0(out.dir.transforms, "/", metadata$Subdirectory[i], "/", out.prefix, "Step2_"))

    tx.img.step2 = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image,
                                       transformlist = step2$fwdtransforms,
                                       interpolator = "linear")
    antsImageWrite(tx.img.step2, paste0(out.dir.txVol, "/", metadata$Subdirectory[i], "/",  out.prefix, "Step2_Warp2Tx_fwdWarpedImg.nrrd"))
    
    
    antsApplyTransforms(fixed= ref.image, 
                        moving = mov.image, 
                        transformlist = step2$fwdtransforms, 
                        compose = paste0(out.dir.transforms, "/", metadata$Subdirectory[i], "/", out.prefix, "combinedFwdWarps-"))
    
    # Inverse labels
    print("Saving inverse-transformed labels")
    inverse.labels1 = antsApplyTransforms(fixed = mov.image, 
                                          moving = ref.label,
                                          transformlist = step1$invtransforms,
                                          whichtoinvert = c(T,F), 
                                          interpolator = "genericLabel")
    antsImageWrite(inverse.labels1, paste0(out.dir.txLab, "/", metadata$Subdirectory[i], "/",out.prefix, "step1-inverse-label.nii.gz"))
    
    
    inverse.labels2 = antsApplyTransforms(fixed = mov.image, 
                                         moving = ref.label,
                                         transformlist = c(step1$invtransforms[1],
                                                           step1$invtransforms[2],
                                                           step2$invtransforms[3]),
                                         whichtoinvert = c(T,F,F), 
                                         interpolator = "genericLabel")
    antsImageWrite(inverse.labels2, paste0(out.dir.txLab, "/", metadata$Subdirectory[i], "/",out.prefix, "step2-inverse-label.nii.gz"))
    
    remove(mov.image, mov.label, mov.mask,
           step1, step2,
           tx.img.step1Linear, tx.img.step1, tx.img.step2,
           inverse.labels1, inverse.labels2)
    gc()
    
    print(paste0("Completed registration ", i, "/", nrow(metadata), ": ", ScanID))
    print(Sys.time()-start)
  }