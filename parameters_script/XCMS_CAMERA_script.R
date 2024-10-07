# XCMS processing and CAMERA annotation 
library(xcms)
library(CAMERA)
library(mzR)
library(stringr)
library(MSnbase)
library("BiocParallel")



run_XCMS_CAMERA<-function(files_dir,
                    bw=5,peakwidth=c(5, 25), ppm=10, snthresh = 100, mzdiff = 0.01,
                    prefilter = c(2, 1000), mzwid = 0.015,noise = 100,
                    perfwhm = 0.6,mzabs = 0.01,cor_eic_th = 0.8,
                    sample_name=sample_name,
                    sample_group =sample_group,
                    return_output="featureTable", # return CAMERA featureTable
                    polarity="positive"){
  
  
  # data import
  dat.dirs<-list.dirs(path=files_dir,recursive = T, full.names = T)
  

  MS_files<-list.files(path=dat.dirs,
                             recursive=F,
                             pattern="mzXML",
                             full.names=T)
  
  if(!is.null(sample_group)) {
    sample_group = sample_group
  }
  else{sample_group =rep(1:length(MS_files))}
  if(!is.null(sample_name)) {
    sample_name = sample_name
  }
  else{
    sample_name =rep(1:length(MS_files))
    }
  
  pd <- data.frame(sample_name =sample_name,
                   sample_group = sample_group,
                   stringsAsFactors = FALSE)
  
  
  raw_data <- readMSData(files=MS_files, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")
  
  ## peak detection
  cwp<-CentWaveParam(peakwidth=peakwidth, noise=noise,
                     prefilter=prefilter)
  xdata <- findChromPeaks(raw_data, param = cwp)
  
  xdata_MS1 <- filterMsLevel(xdata, msLevel = 1L)
  
  ## Alignment
  xdata_MS1<- adjustRtime(xdata_MS1, param = ObiwarpParam())
  
  ## Perform the correspondence
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                          minFraction = 0.5, bw = bw)
  
  xdata_MS1 <- groupChromPeaks(xdata_MS1, param = pdp)
  xdata_MS1 <- fillChromPeaks(xdata_MS1, param = ChromPeakAreaParam())


  ## Annotation with CAMERA

  xset <- as(xdata_MS1, 'xcmsSet')
  ## Create an xsAnnotate object
  an <- xsAnnotate(xset)
  ## Group based on RT
  anF <- groupFWHM(an, perfwhm = perfwhm)
  ## Annotate isotopes
  anI <- findIsotopes(anF, mzabs =mzabs)
  ## Verify grouping
  anIC <- groupCorr(anI,  cor_eic_th= cor_eic_th)
  ## Annotate adducts
  anFA <- findAdducts(anIC, polarity=polarity)
  featureTable <- getPeaklist(anFA)
  featureTable<-featureTable[order(featureTable$rt),]
  rownames(featureTable)<-paste0("FT",rownames(featureTable))
  return(featureTable)
}


## Test: read in with demo mzXML files
run_XCMS_CAMERA_test<-run_XCMS_CAMERA(files_dir = "C:/Users/Downloads/Si11_mzXML/Si11_mzXML_files",
                                      sample_name=c(paste("sample", 1:7)),sample_group=c(paste("sample", 1:7)))

