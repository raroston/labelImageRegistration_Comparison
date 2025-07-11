# Tensory-Based Morphometry (TBM)
# Rachel Roston, Ph.D.

library(ANTsR)

# intensity-only TBM
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

# TBM
  ## setup output directory
  dir.out.tbm = paste0("./Results/", experiment, "/TBM-", txPattern)
  if(!dir.exists(dir.out.tbm)) dir.create(dir.out.tbm)
  
 ## calculate Jacobians
  warps = lapply(paste0(tx.dir, "/", tx.paths), antsImageRead)
  jacobians = lapply(warps, createJacobianDeterminantImage, domainImg = ref.img, geom = T)
  
  ## analysis
    #set the regression
    print("Starting statistical analysis")
    j_mat = imageListToMatrix( jacobians, ref.mask )
    j_mdl = lm( j_mat ~ genotypes)    # raw volume data
    j_bmdl = bigLMStats( j_mdl , 1.e-5 )
    print( min( p.adjust( j_bmdl$beta.pval, 'none' ) , na.rm=T ) )
    
    #do multiple corrections
    #P defines the significance level we want to evaluate things after FDR correction
    volstats = j_bmdl$beta.pval
    
    pAdjMethod = "fdr"
    correct.ps= p.adjust(volstats, pAdjMethod)
    pImg = makeImage(ref.mask, correct.ps)
    
    P = 0.05
    pMask = thresholdImage(pImg, 10^-16, P, 1, 0)
    betaImg = makeImage( ref.mask, j_bmdl$beta ) * pMask
    
    print("Writing images")
    antsImageWrite(betaImg, paste0(dir.out.tbm, "/", experiment, "_", pAdjMethod, "_betaImg.nii.gz"))
}
    
# label-intensity TBM   
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

# TBM
  ## setup output directory
  dir.out.tbm = paste0("./Results/", experiment, "/TBM-", txPattern)
  if(!dir.exists(dir.out.tbm)) dir.create(dir.out.tbm)
  
  ## calculate Jacobians
  warps = lapply(paste0(tx.dir, "/", tx.paths), antsImageRead)
  jacobians = lapply(warps, createJacobianDeterminantImage, domainImg = ref.img, geom = T)
  
  ## analysis
    #set the regression
    print("Starting statistical analysis")
    j_mat = imageListToMatrix( jacobians, ref.mask )
    j_mdl = lm( j_mat ~ genotypes)    # raw volume data
    j_bmdl = bigLMStats( j_mdl , 1.e-5 )
    print( min( p.adjust( j_bmdl$beta.pval, 'none' ) , na.rm=T ) )
    
    #do multiple corrections
    #P defines the significance level we want to evaluate things after FDR correction
    volstats = j_bmdl$beta.pval
    
    pAdjMethod = "fdr"
    correct.ps= p.adjust(volstats, pAdjMethod)
    pImg = makeImage(ref.mask, correct.ps)
    
    P = 0.05
    pMask = thresholdImage(pImg, 10^-16, P, 1, 0)
    betaImg = makeImage( ref.mask, j_bmdl$beta ) * pMask
    
    print("Writing images")
    antsImageWrite(betaImg, paste0(dir.out.tbm, "/", experiment, "_", pAdjMethod, "_betaImg.nii.gz"))
} 