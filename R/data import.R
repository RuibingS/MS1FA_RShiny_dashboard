
## read feature table csv file

featureTable.import.fun <- function(file.dir) {
  
  featureTable <- tryCatch({
    
    read.csv(file.dir, check.names = TRUE)
  }, error = function(e) {
    
    tryCatch({
      read.csv2(file.dir, check.names = TRUE)
    }, error = function(e2) {
     
      stop(paste("Error: Unable to read the CSV file. Details:", e2$message))
    })
  })
  
 
  if (nrow(featureTable) == 0 || ncol(featureTable) == 1) {
    stop("Error: Unable to read the CSV file correctly. Please make sure it uses a valid format CSV.")
    
  }

  
  #featureTable <- read.csv(file.dir,check.names = TRUE)
  FT_colnames <- colnames(featureTable)

  # change column names
  # check if "feature_name", change "mzmed" to "mz", "rtmed" to "rt"
  tryCatch({
    # Check if column names are empty or do not include "mz" and "rt"
    if (any(colnames(featureTable) == "")) {
      stop("Column names cannot be empty.")
    }
    if (!any(c("mz", "rt","mzmed","rtmed","row m/z","row retention time","row.m.z","row.retention.time") %in% colnames(featureTable))) {
      stop("The 'mz' and 'rt' columns are required.")
    }
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
  })
 
  # xcms online output: feature_name starting with M1111T18
  if("feature_name"%in% FT_colnames){
    if (grepl("M", featureTable[1, "feature_name"])){
      # remove the feature_name 
      featureTable_sub<-featureTable[,-which(colnames(featureTable)%in%c("feature_name"))]
      
      if (is.unsorted(featureTable_sub$rt)) {
        featureTable_order <- featureTable_sub[order(featureTable_sub$rt), ]
        featureTable_order$feature_name <- sprintf("FT%04d", seq(1, nrow(featureTable_order)))
        
        
        column_to_move <- "feature_name"
        desired_index <- 1  
        featureTable_order_df <- featureTable_order %>%
          relocate(!!rlang::sym(column_to_move), .before = 1)
        
        featureTable_order_df$mz <- as.numeric(featureTable_order_df$mz)
        featureTable_order_df$rt <- as.numeric(featureTable_order_df$rt)
        
        return(featureTable_order_df)
    }
      else{
      featureTable_sub$feature_name <- sprintf("FT%04d", seq(1, nrow(featureTable_sub)))
      column_to_move <- "feature_name"
      desired_index <- 1  
      featureTable_order_df <- featureTable_sub %>%
        relocate(!!rlang::sym(column_to_move), .before = 1)
      
      featureTable_order_df$mz <- as.numeric(featureTable_order_df$mz)
      featureTable_order_df$rt <- as.numeric(featureTable_order_df$rt)
      return(featureTable_order_df)
      }
    }
  }
  
  # xcms output:
  if (any(!grepl("feature_name", FT_colnames) & grepl("FT", featureTable[1,1])) & !"row.ID"%in%FT_colnames) 
    {
   
    colnames(featureTable)[1] <- "feature_name"
    
    if ("mzmed" %in% FT_colnames) {
      colnames(featureTable)[colnames(featureTable) == "mzmed"] <- "mz"
    }
    if ("rtmed" %in% FT_colnames) {
      colnames(featureTable)[colnames(featureTable) == "rtmed"] <- "rt"
    }
    
    featureTable$mz <- as.numeric(featureTable$mz)
    featureTable$rt <- round(as.numeric(featureTable$rt), 3)
  
    if (is.unsorted(featureTable$rt)) {
    
      featureTable$feature_name[]<-""
      featureTable_order <- featureTable[order(featureTable$rt), ]
      featureTable_order$feature_name <- sprintf("FT%04d", seq(1, nrow(featureTable_order)))
      
      return(featureTable_order)
    }
    else{
     
      return(featureTable)
    }
  } else if (!"feature_name" %in%FT_colnames & (!"id" %in%FT_colnames)  & !"row.ID"%in%FT_colnames )
    {

    if (is.unsorted(featureTable$rt) ) {
      # sort by RT
      featureTable_order <- featureTable[order(featureTable$rt), ]
      # add a new column 'feature_name' with formatted values
      featureTable_order$feature_name <- sprintf("FT%04d", seq(1, nrow(featureTable_order)))
      
      # move to the first index
      column_to_move <- "feature_name"
      desired_index <- 1
      featureTable_order_df <- featureTable_order %>%
        relocate(!!rlang::sym(column_to_move), .before = 1)
      featureTable_order_df$mz <- as.numeric(featureTable_order_df$mz)
      featureTable_order_df$rt <- as.numeric(featureTable_order_df$rt)
      
      return(featureTable_order_df)
    }
    else{
     
      featureTable$feature_name <- sprintf("FT%04d", seq(1, nrow(featureTable)))
      column_to_move <- "feature_name"
      desired_index <- 1  
      featureTable_order_df <- featureTable %>%
        relocate(!!rlang::sym(column_to_move), .before = 1)

      featureTable_order_df$mz <- as.numeric(featureTable_order_df$mz)
      featureTable_order_df$rt <- as.numeric(featureTable_order_df$rt)
      return(featureTable_order_df)
    }
    
  }  
  # MZmine feature table
  if ("row.ID" %in% FT_colnames) { 

      
    colnames(featureTable)[1] <- "feature_name"
    colnames(featureTable)[2] <- "mz"
    colnames(featureTable)[3] <- "rt"
    featureTable$feature_name<-sprintf("FT%04d",featureTable$feature_name)
    featureTable$mz <- as.numeric(featureTable$mz)
    featureTable$rt <- round(as.numeric(featureTable$rt) * 60, 3)

    return(featureTable)
  
    
  }
  
}

# Check if RT is in seconds if not convert to seconds
check_rt_column <- function(FT_data) {
  rt_message <- "" 
  
  if ("rt" %in% colnames(FT_data)) {
    max_rt <- max(FT_data$rt, na.rm = TRUE)
    if (max_rt <= 60) {
     
      FT_data$rt <- FT_data$rt * 60
      rt_message <- "rt in seconds needed, converted from minutes to seconds."
    }
  }
  
  return(list(data = FT_data, message = rt_message))
}


######################################################################################
# read MS2 .mgf file from MZmine output 

MS2_mgf_import_fun<-function(mgf_file){
  MS2_mgf_data<-read.delim(mgf_file,header=FALSE)
  MS2_mgf_list<-split(MS2_mgf_data[1:nrow(MS2_mgf_data), ], cumsum(1:nrow(MS2_mgf_data)%in%c(grep("^BEGIN IONS", MS2_mgf_data[1:nrow(MS2_mgf_data), ]))))

  
  parseStringToMS2<-function(data){
    if(length(data)> 0 ){
      #set MS2Df
      MS2Df<-NA
      if(sum(grepl('^BEGIN IONS', data)==T)>0){
        MS2Df<-data.frame(matrix(as.numeric(unlist(strsplit(data[(which(grepl('^MSLEVEL=2', data) == TRUE)+1):(length(data)-1)], " ", fixed=FALSE))),ncol=2, byrow=TRUE))
        colnames(MS2Df)<-c("mz","intensity")
      
        }
      
    #set FEATURE_ID
    feature_name<-ifelse(any(grepl("^FEATURE_ID", data)==TRUE),sprintf("FT%04d",as.numeric(unlist(stringr::str_split(data[grep("FEATURE_ID",data)], "=",2))[2])),NA_character_)
    
    #set PEPMASS
    PI_mass<-as.numeric(ifelse(any(grepl("^PEPMASS", data)==TRUE),unlist(stringr::str_split(data[grep("PEPMASS",data)], "=",2))[2],NA_real_))
    #set RT
    rt<-as.numeric(ifelse(any(grepl("^RTINSECONDS", data)==TRUE),unlist(stringr::str_split(data[grep("RTINSECONDS",data)], "=",2))[2],NA_real_))
    
    #set scan
    scans<-as.numeric(ifelse(any(grepl("^SCANS", data)==TRUE),unlist(stringr::str_split(data[grep("SCANS",data)], "=",2))[2],NA_real_))
    
    
    list(
      feature_name=feature_name,
      PI_mass=PI_mass,
      scans=scans,
      rt=rt,
      MS2Df=MS2Df
    )
  }else{
    NULL
  }
  }
    
  libOutputFun<-function(libDataInput){
    
    libDataOutput<-lapply(libDataInput, function(x) parseStringToMS2(x))
    ## remove NULL
    libDataOutput[sapply(libDataOutput, is.null)==FALSE]
  }
  
  ##Pipeline: Import data and parse it into a List of compounds
  output<-libOutputFun(MS2_mgf_list)

  return(output)
  
}



#######################################################################################
## function import metabolite data- .csv

metabolie.data.import.fun<-function(metabo_data.path,IonPolarity){

  metabo_data<-read.csv(metabo_data.path)

  #remove columns where all vales are NA
  metabo_data <- metabo_data[,colSums(is.na(metabo_data))<nrow(metabo_data)]

  # calulate exact mass from formula
  exactMass<-vector()
  for(i in 1:nrow(metabo_data)){
    exactMass[i]<- isopattern(isotopes,
                              metabo_data[i,"Formula"],
                              threshold=0.1,
                              plotit=F,
                              charge=FALSE,
                              emass=0.00054858,
                              algo=1)[[1]][1]

  }
  #calculate primary Adducts "[M+H]+","[M+Na]+"
  if(IonPolarity =="pos"){
    metabo_data_add<-MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M+H]+","[M+Na]+"))
  }

  #calculate primary Adducts [M-H]-
  if(IonPolarity =="neg"){
      metabo_data_add<-MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M-H]-"))#[,1]
  }


  return(cbind(metabo_data,exactMass,metabo_data_add))
}
#metabo_data<-metabolie.data.import.fun("~/Downloads/MS1FA_RShiny_dashboard/Data/metabolite_data/Si11_target_list.csv",IonPolarity="pos")
#metabo_data<-metabolie.data.import.fun("~/Downloads/MS1FA_RShiny_dashboard/Data/metabolite_data/target_list_01_26_2024.csv",IonPolarity="pos")

#######################################################################################
## read .library file
######################################

read_library_Fun<-function(lib_dir,IonPolarity ){
  
  ## if multiple libraries, then concatenate them as lists
  read_lib_fun<-function(lib_dir){
    lib.data<-read.delim(lib_dir,header=FALSE)
    return(lib.data)
  }
  
  
  parseStringToCompound<-function(data){
    if(length(data)> 0 ){
      #set spectraDf
      spectraDf<-NA
      if(sum(grepl('^Num Peaks', data)==T)>0){
        spectraDf<-data.frame(matrix(as.numeric(unlist(strsplit(data[(which(grepl('^Num Peaks', data) == TRUE)+1):length(data)], " ", fixed=FALSE))),ncol=2, byrow=TRUE))
        colnames(spectraDf)<-c("mz","rel_intensity")
      }
      
      #set formula
      Formula<-ifelse(any(grepl("^Formula", data)==TRUE),unlist(stringr::str_split(data[grep("Formula",data)], ": ",2))[2],NA_character_)
      
      #set compound name
      Name<-unlist(stringr::str_split(data[1], ": ",2))[2]
      
      #set ionPolarity
      ionPolarity<-unlist(stringr::str_split(data[grep("IonPolarity",data)], ": ",2))[2]
      
      #set msLevel
      msLevel<-if (sum(grepl("^MSMS", data))==0) {paste("MS1")}
      else{ifelse(sum(grepl("^MSMS: 1", data)==T)>0, paste("MS1"),paste("MS2"))}
      
      #set rt
      RT<-ifelse(any(grepl("RetTime",data)), as.numeric(unlist(stringr::str_split(data[[grep("RetTime",data)]], ": ",2))[2]),NA_real_)
      
      if(!is.na(Formula) ){
        exactMass<-ifelse(!is.na(Formula),
                          isopattern(isotopes,
                                     Formula,
                                     threshold=0.1,
                                     plotit=F,
                                     charge=FALSE,
                                     emass=0.00054858,
                                     algo=1)[[1]][1],
                          NaN)
        #calculate primary adducts
        if(IonPolarity =="pos"){
          primaryAdducts<-MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M+H]+","[M+Na]+","[M+K]+","[M+NH4]+"))
          primaryAdducts<-cbind.data.frame(Name,exactMass,Formula,RT,primaryAdducts) #,check.names=FALSE
          
        }
        if(IonPolarity =="neg"){
          primaryAdducts<-MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M-H]-"))
          primaryAdducts<-cbind.data.frame(Name,exactMass,Formula,RT,primaryAdducts) #,check.names=FALSE
          
        }
        
        
        list(
          name=Name,
          formula=formula,
          ionPolarity=ionPolarity,
          msLevel=msLevel,
          rt=rt,
          spectraDf=spectraDf,
          exactMass=exactMass,
          primaryAdducts=primaryAdducts
        )
      }else{
        NULL
      }
    }else{
      NULL
    }
  }
  
  ## --Step 1--
  ## read the '.library' data and returns a list of vectors
  readDataFun<-function(lib_dir){
    libData<-read_lib_fun(lib_dir)
   
    splitData<-split(libData[1:nrow(libData), ], cumsum(1:nrow(libData)%in%c(grep("^Name", libData[1:nrow(libData), ]))))
    splitData
  }
  
  ## --Step 2--
  ## loop over dataset and parse the data into a list of objects
  libOutputFun<-function(libDataInput){
    
    libDataOutput<-lapply(libDataInput, function(x) parseStringToCompound(x))

    libDataOutput[sapply(libDataOutput, is.null)==FALSE]
  }
  
  
  ##Pipeline: Import data and parse it into a List of compounds
  output<-libOutputFun(readDataFun(lib_dir))
  output_names<-sapply(output, function(x) x$name) ## remove duplicated names
  output<-output[which(!duplicated(output_names))]
  
  output.df<-do.call(rbind,lapply(output, function(x) x$primaryAdducts))
  return(output.df)
}


# test
#library_teest<-read_library_Fun(lib_dir = "~/Downloads/MS1FA_RShiny_dashboard/Data/metabolite_data/PA14.library",
#                                IonPolarity = "neg")
#library_teest<-read_library_Fun(lib_dir = "~/Downloads/MS1FA_RShiny_dashboard/Data/metabolite_data/PA14.library",
#                                IonPolarity = "pos")
#head(library_teest)
