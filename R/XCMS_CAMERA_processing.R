# # Si11 spike in PA14 samples
# 
# library(xcms)
# 
# library(RColorBrewer)
# 
# library(pheatmap)
# 
# 

# files_dir_spikein <- "~/Downloads/MS/Si16_pool_mixture"
# 
# MS_files <- list.files(files_dir_spikein,
#                        recursive = FALSE,
#                        pattern = "\\.mzXML$|\\.mzML$",
#                        full.names = TRUE)
# 
# 
# 
# pd <- data.frame(sample_name = c("StM16_10uM_A","StM16_10uM_B","StM16_10uM_C",
#                                  "StM16_4uM_A","StM16_4uM_B","StM16_4uM_C",
#                                  "StM16_8uM_A", "StM16_8uM_B" ,"StM16_8uM_C"),
#                  sample_group = c(rep("StM16_10uM",3),rep("StM16_4uM",3),
#                                   rep("StM16_8uM",3)
#                  ),
#                  stringsAsFactors = FALSE)


#
# # Load raw MS data
# raw_data <- MSnbase::readMSData(files = MS_files,
#                                 pdata = new("NAnnotatedDataFrame", pd),
#                                 mode = "onDisk")

#  # Load raw MS data
# raw_data <- MSnbase::readMSData(files = MS_files,
#                                  pdata = new("NAnnotatedDataFrame", pd),
#                                  mode = "onDisk")
#  
# 
# cwp <- CentWaveParam(
#    peakwidth = c(5, 12),      
#    ppm = 12,                 
#    snthresh = 15,            
#    mzdiff = 0.02,            
#    prefilter = c(4, 1000),   
#    noise = 1000               
#  )
#  
# xdata <- findChromPeaks(raw_data, param = cwp)
#  
# 
# xdata_MS1 <- filterMsLevel(xdata, msLevel = 1L)
#  
# 
# xdata_MS1 <- adjustRtime(xdata_MS1, param = ObiwarpParam())
# 
# pdp <- PeakDensityParam(
#    sampleGroups = xdata$sample_group,
#    minFraction = 0.7,  
#    bw = 2.0            
#  )
#  
# xdata_MS1 <- groupChromPeaks(xdata_MS1, param = pdp)
#  
#  
# xdata_MS1 <- fillChromPeaks(xdata_MS1, param = ChromPeakAreaParam())
#  
# 
# fea_def <- featureDefinitions(xdata_MS1)
# fea_int <- featureValues(xdata_MS1, value = "into")
# fea_output <- cbind(fea_def, fea_int)
#  
# 
# dim(fea_output)  
# 
# write.csv(fea_output[,-which(colnames(fea_output)%in%"peakidx")], "Spike_in_Si16_XCMS.csv")
# 
# library(CAMERA)
# 
# xset <- as(xdata_MS1, 'xcmsSet')
# xsa <- xsAnnotate(xset)
# 
############################################################################################################
#  CAMERA annotation                                                                                       #
############################################################################################################ 
# anF <- groupFWHM(xsa, perfwhm = 0.6)
# 
# anI <- findIsotopes(anF, mzabs = 0.01)
# 
# anIC <- groupCorr(anI, cor_eic_th = 0.75)
# 
# anFA <- findAdducts(anIC, polarity="positive")
# anFA
# 
# peaklist <- getPeaklist(anFA)
# 
# head(peaklist)
# 
# #write.csv(peaklist, file='CAMERA_annotated_FT_31012025_MS1.csv')
# 
# 
# 
# write.csv(peaklist, file='CAMERA_annotated_FT_01022025_MS1_2500.csv')
# 
# 







