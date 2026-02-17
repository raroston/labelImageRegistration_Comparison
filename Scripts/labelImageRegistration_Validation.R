library(ANTsR)
Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = 24)

out.dir = "./Results/validation_output/"
dir.create(out.dir)

out.table = "./Results/label_overlaps.csv" 
if(!file.exists(out.table)) file.create(out.table)

metadata = read.csv(file="./Data/scan_genotypes.csv")

ko = which(metadata$Genotype == "MUT/MUT")
wt = which(metadata$Genotype == "WT/WT")

metadata = metadata[c(wt, ko), ]

ref.image = antsImageRead("./Data/Template/maskedtemplate0__lowRes.nii.gz")
ref.label = antsImageRead("./Data/Template/maskedtemplate0__lowRes-15-label.nii.gz")
ref.mask  = antsImageRead("./Data/Template/maskedtemplate0__lowRes-wholebody-label.nii.gz")

n=nrow(metadata)
results=NULL

label.order = c(4,5,3,10,15,7,8,11,2,1,6,12,13,14,9)

for (label in label.order) {
  for (i in 1:n) {
    if(!file.exists(paste0(out.dir, "/", label, "/", metadata$ScanID[i], "-label.nii.gz"))){
      print(paste0("Starting specimen: ", i, " of ", n))
      intensityOnly.label = antsImageRead(dir(patt=paste0(metadata$ScanID[i], "__intensityOnly_step2"), path = './Results/intensityOnly/invTx_Labels/', recursive = TRUE, full.names = TRUE))
      mov.image = antsImageRead(dir(patt=metadata$ScanID[i], path='./Data/Volumes/', full.names = TRUE) )
      mov.label = antsImageRead(dir(patt=metadata$ScanID[i], path='./Data/Labels/', full.names = TRUE) )
      mov.mask = antsImageRead(dir(patt=metadata$ScanID[i], path='./Data/Masks/', full.names = TRUE) )
    
      # begin intensity registration leave one out 
    
      print(paste0("starting label ", label))
      ref.label.tmp = antsImageClone( ref.label)
      ref.label.tmp[ref.label==label] = 0 #remove the label being held.
      mov.label.tmp = antsImageClone( mov.label)
      mov.label.tmp[mov.label==label] = 0 #remove the label being held.
      
      # Step1
      start.time = Sys.time()
      print(paste0("Starting step1, ", Sys.time()))
      
      step1 = labelImageRegistration(fixedLabelImages = ref.label.tmp,
                                     movingLabelImages = mov.label.tmp,
                                     fixedIntensityImages = ref.image,
                                     movingIntensityImages = mov.image,
                                     fixedMask = ref.mask,
                                     movingMask = mov.mask,
                                     initialTransforms = "similarity",
                                     typeOfDeformableTransform = "antsRegistrationSyN[so]",
                                     outputPrefix = paste0(out.dir, "tmp-step1-"))
      print(paste0("Finished step1: ", Sys.time()-start.time))
      
      #Step2
      start.time = Sys.time()
      paste0("Starting step2, ", Sys.time())
      step2 = antsRegistration(fixed = ref.image,
                               moving = mov.image,
                               typeofTransform = "antsRegistrationSyN[so]",
                               mask = ref.mask,
                               movingMask = mov.mask,
                               initialTransform = step1$fwdtransforms, 
                               outprefix = paste0(out.dir, "tmp-step2-"))
      print(Sys.time()-start.time)
      
      # Inverse labels
      ref.labels.in.subject.label = antsApplyTransforms(fixed = mov.image, 
                                            moving = ref.label,
                                            transformlist = c(step1$invtransforms[1],
                                                              step1$invtransforms[2],
                                                              step2$invtransforms[3]),
                                            whichtoinvert = c(T,F,F), 
                                            interpolator = "genericLabel")
      if(!dir.exists(paste0(out.dir, label))) dir.create(paste0(out.dir, label))
      antsImageWrite(ref.labels.in.subject.label, paste0(out.dir, label,"/", metadata$ScanID[i], "-label.nii.gz"))
      
      # make new images with just the held label in them
      manual.label = antsImageClone(mov.label)
      manual.label[mov.label!=label] = 0
      manual.label[mov.label==label] = 1
      
      int.label = antsImageClone(intensityOnly.label)
      int.label[intensityOnly.label != label] = 0
      int.label[intensityOnly.label == label] = 1
      
      label.label = antsImageClone(ref.labels.in.subject.label)
      label.label[ref.labels.in.subject.label != label] = 0
      label.label[ref.labels.in.subject.label == label] = 1
      
      overlap.int = labelOverlapMeasures(int.label , manual.label)[1,]
      overlap.label = labelOverlapMeasures(label.label , manual.label)[1,]
      
      results = cbind(ScanID = metadata$ScanID[i], Label.Omit = label, overlap.int, overlap.label)
      if(length(readLines(out.table)) < 1) {
        write.table(results, 
                    file = out.table, 
                    append = F, sep = ",",
                    row.names = F,
                    col.names = T)
      } else {
        write.table(results, 
                    file = out.table, 
                    append = T, sep = ",",
                    row.names = F,
                    col.names = F)
      }
      
      print(Sys.time()-start.time)
      
      remove(mov.image, mov.label, ref.label.tmp, mov.label.tmp, step1, step2, ref.labels.in.subject.label, manual.label, int.label, label.label, overlap.label, overlap.int)
      gc()
      }
  }
}
