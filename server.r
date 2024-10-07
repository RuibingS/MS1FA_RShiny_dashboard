library(BiocManager)
library(igraph)
library(stringr)
library(readxl)
library(purrr)
library(readr)
library(MetaboCoreUtils) # for function: mass2mz
library(enviPat) # for function: isopattern
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(viridisLite)
library(ggplot2)
library(roxygen2)
library(rlang)
library(RcppArmadillo)
library(pryr)
library(webshot)
library(htmlwidgets)
library(profvis)
library(magick)
library(bs4Dash)
library(MSnbase)
library(chromote)
library(mzR)
library(here)


source(here::here("R","required_packages.R"))
data(isotopes)

source(here::here("R","data import.R"))

source(here::here("R","helper functions.R"))



Rcpp::sourceCpp(here::here("src",'findMatch.cpp'))
Rcpp::sourceCpp(here::here("src",'MS2_mzMatch.cpp'))
Rcpp::sourceCpp(here::here("src",'MS2_precursor_match.cpp'))
Rcpp::sourceCpp(here::here("src",'IsfAnno.cpp'))
Rcpp::sourceCpp(here::here("src",'Check_iso_charges_function.cpp'))
Rcpp::sourceCpp(here::here("src",'GroupIndex.cpp'))
Rcpp::sourceCpp(here::here("src",'pairwise_CorFunction.cpp'))
Rcpp::sourceCpp(here::here("src",'cor_index.cpp'))
Rcpp::sourceCpp(here::here("src","adducts_match.cpp"))
Rcpp::sourceCpp(here::here("src","Isfanno_mzmine.cpp"))
Rcpp::sourceCpp(here::here("src","left_join_rcpp.cpp"))


options(shiny.maxRequestSize=50*1024^2)


default_adduct_file_path<-here::here("Data", "Adducts","adducts.csv")
default_NL_file_path <- here::here("Data", "Neutral loss", "neutral loss database.csv")

# demo files
# feature table and MS2- PA14 
demoPA14_FT_path<-here::here("Data", "feature_table","PA14_featureTable_XCMS_CAMERA.csv") 
demo_PA14_MS2_path<-here::here("Data", "MS2_files","MS2_antibiotic","p_hm188_pool_25745.mzXML") # pool sample

# feature table and MS2- Si11 
demoSi11_FT_path<-here::here("Data", "feature_table","si11_MZmine3_iimn_gnps_quant.csv") # MZmine FT table
demo_Si11_MS2_path<-here::here("Data", "MS2_files","Si11_MS2_MZmine3_mgf","si11_MZmine3_iimn_gnps.mgf") # MZmine mgf file

# target list
demo_PA14Metabo_path<-here::here("Data", "metabolite_data","target_list_V1_2024.csv") 
demo_Si11Metabo_path<-here::here("Data", "metabolite_data","Si11_target_list.csv")



server <- function(input, output,session) {
 

  output$memoryUsage <- renderText({
    invalidateLater(5000, session)  # Invalidate this reactive expression every 5 seconds
    paste("Memory used (MB):", round(pryr::mem_used() / 10^6, 2))
  })

  observe({
    
    minValue <- min(input$minValue, input$maxValue)
    maxValue <- max(input$minValue, input$maxValue)
    
    # Update the sliderInput
    updateSliderInput(session, "RTrange",
                      value = c(minValue, maxValue))
  })
  # Demo feature table - csv file
  
  

  

  selFilePath <- reactiveVal(NULL)
  fileSource <- reactiveVal("none")  
  uploadFTFile <- reactiveVal("No file selected")
  
  observeEvent(input$uploadDemoTable, {
    if (!is.null(input$demoFTFile)) {
      if (input$demoFTFile == "PA14example1") {
        selFilePath(demoPA14_FT_path)
        fileSource("demo")
        uploadFTFile("Demo file selected")
      } else if (input$demoFTFile == "Si11example2") {
        selFilePath(demoSi11_FT_path)
        fileSource("demo")
        uploadFTFile("Demo file selected")
      }
     
    }
  })
  

  observeEvent(input$FT_file, {
    selFilePath(input$FT_file$datapath)
    fileSource("user")
    uploadFTFile(input$FT_file$name)
   
  })
  
  
  
  # select default feature table to run
  output$FTfile <- renderUI({
    fileInput(
      inputId = "FT_file", 
      label = "Choose a feature table file (.csv)",
      placeholder = uploadFTFile(),
      multiple = FALSE
    )
  })
  

  output$featureTable <- renderDataTable({
    req(input$showTable)
    
    FT_table<-get.df_input()
    FT_table_check<-check_rt_column(FT_table)$data
    datatable(FT_table_check)
  })

  # Render the demo feature table based on the checkbox input
  output$demoFeatureTable <- renderDataTable({
    req(input$viewDemoTable)
    
    
    if (input$viewDemoTable) {
      demoFTFile <- input$demoFTFile
      if (is.null(demoFTFile)) {
        demoFTFile <- "PA14example1" 
      }
      
      path <- if (demoFTFile == "PA14example1") demoPA14_FT_path else demoSi11_FT_path
      
      data <- read.csv(path)
      datatable(data)
    }
  })
  
  # Download handlers for example csv file - PA14
  output$download1 <- downloadHandler(
    filename = function() {
      "Demo_PA14_featureTable_XCMS.csv"
    },
    content = function(file) {
      file.copy(demoPA14_FT_path, file)
    }
  )
  # Download handlers for example csv file - Si11
  output$download2 <- downloadHandler(
    filename = function() {
      "Demo_Si11_featureTable_MZmine.csv"
    },
    content = function(file) {
      file.copy(demoSi11_FT_path, file)
    }
  )

  
  # feature table - csv file
  get.df_input<- reactive({
    req(selFilePath())
    path<-selFilePath()
    tryCatch({
      FT <- featureTable.import.fun(path)
      
      return(FT)
    }, error = function(e) {

      output$error_message <- renderText({ e$message })
      NULL
    })
    
  })

  
   get.df<- reactive({
    FT_table<-get.df_input()
    FT_table_check<-check_rt_column(FT_table)$data
   
    FT_sub<-FT_table_check %>% #
      dplyr::filter(rt >= input$RTrange[1] & 
                      rt <= input$RTrange[2])
    return(FT_sub)
  })

  
  ########## MS2 file
  selectedMS2FilePath <- reactiveVal(NULL)
  
  get.ms2df <- reactive({
    req(selectedMS2FilePath())
    if (!file.exists(selectedMS2FilePath())) {
      message("Selected MS2 file path does not exist. Skipping MS2 processing.")
      return(NULL)  
    }
    if(any(grepl("mzXML",selectedMS2FilePath())|grepl("mzML",selectedMS2FilePath())))
      {
      validate(need(selectedMS2FilePath(), "Please upload a pooled MS2 file()"))
      MSnbase::readMSData(files = selectedMS2FilePath(),
                          
                           msLevel. = 2, mode = "onDisk")
    

  
      }
    if(any(grepl(".mgf",selectedMS2FilePath())))
      {
      MS2_mgf_import<-MS2_mgf_import_fun(mgf_file=selectedMS2FilePath())
      MS2_mgf_import
     
      }
  })
  
 
  
  fileSource_MS2 <- reactiveVal("none")  
  uploadedMS2FileName <- reactiveVal("No file selected")
  
  observeEvent(input$UploadMS2DemoTable, {
    if (!is.null(input$demoMS2File)) {
      if (input$demoMS2File == "MS2example1") {
        selectedMS2FilePath(demo_PA14_MS2_path)
      } else if (input$demoMS2File == "MS2example2") {
        selectedMS2FilePath(demo_Si11_MS2_path)
      }
      fileSource_MS2("demo")
      uploadedMS2FileName("Demo file selected")

    }
  })
  
  observeEvent(input$MS2_file, {
    
    selectedMS2FilePath(input$MS2_file$datapath)
    fileSource_MS2("user")
    uploadedMS2FileName(input$MS2_file$name)

  })
  
  # select default feature table to run
  output$MS2file <- renderUI({
    fileInput(
      inputId = "MS2_file", 
      "Choose a MS2 file (.mzXML, .mzXL or .mgf)", 
    placeholder = uploadedMS2FileName(),
    multiple = FALSE)
  })
  
  output$demoMS2Table <- renderDataTable({
    req(input$viewMS2DemoTable)
   
    
    if(input$viewMS2DemoTable){
      
      if(any(grepl("mzXML",demo_PA14_MS2_path)))
      { 
        path <- demo_PA14_MS2_path
        MS2data<- MSnbase::fData(MSnbase::readMSData(files = path, msLevel. = 2, mode = "onDisk"))
        return(MS2data) 
      }
      if(any(grepl(".mgf",demo_Si11_MS2_path)))
      {
        path <- demo_Si11_MS2_path
        MS2_mgf_import<-MS2_mgf_import_fun(mgf_file=path)
        return(plyr::ldply(MS2_mgf_import, data.frame))
      }
    }
  })

  
  output$MS2table <- renderDataTable({
    
    req(input$showMS2Table,selectedMS2FilePath())

    path <- selectedMS2FilePath()
    if(any(grepl("mzXML",path)|grepl("mzML",path)))
    {
      MS2data<- MSnbase::fData(MSnbase::readMSData(files = path, msLevel. = 2, mode = "onDisk"))
      return(MS2data) 

    }
    if(any(grepl(".mgf",path)))
    {
      MS2_mgf_import<-MS2_mgf_import_fun(mgf_file=path)
      return(plyr::ldply(MS2_mgf_import, data.frame))
      
    }
  })
  
  # Download handlers for PA14 MS2 mzXML example file
  output$download3 <- downloadHandler(
    filename = function() {
      "Demo_PA14_MS2.mzXML"
    },
    content = function(file) {
      file.copy(demo_PA14_MS2_path, file)
    }
  )
  # Download handlers for Si11 MS2 mgf example file
  output$download4 <- downloadHandler(
    filename = function() {
      "Demo_Si11_MS2.mgf"
    },
    content = function(file) {
      file.copy(demo_Si11_MS2_path, file)
    }
  )
  
  

  ###################################################################
  

  selTargetFilePath <- reactiveVal(NULL)
  fileSource_metabo <- reactiveVal("none") 
  uploadMetaboFileName <- reactiveVal("No file selected")
    

  # .csv file
  get.metabo<-reactive({

    if(!is.null(selTargetFilePath())){
      # the uploaded file path
      target_path <- selTargetFilePath() 
      
      if (stringr::str_ends(target_path, "csv")) {
        return(metabolie.data.import.fun(target_path, IonPolarity = input$IonPolarity))
      } else if (stringr::str_ends(target_path, "library")) {
        return(read_library_Fun(lib_dir = target_path, IonPolarity = input$IonPolarity))
      }
      
    }
    else{
      return(NULL)
      print(paste0("is.null(selTargetFilePath()",is.null(selTargetFilePath())))
    }
    
  })
  

    observeEvent(input$UploadDemoMetabo, {
      
      if (!is.null(input$demoMetaboFile)) {
        if (input$demoMetaboFile == "PA14_target_example1") {
          selTargetFilePath(demo_PA14Metabo_path)
        } else if (input$demoMetaboFile == "Si11_target_example2") {
          selTargetFilePath(demo_Si11Metabo_path)
        }
        fileSource_metabo("demo")
        uploadMetaboFileName("Demo file selected")
      }
     
    })
    
    observeEvent(input$metabo_file, {
      selTargetFilePath(input$metabo_file$datapath)
      fileSource_metabo("user")
      uploadMetaboFileName(input$metabo_file$name)
    })
    
    # Output for the file input
    output$metabolites_file <- renderUI({
      fileInput(
        inputId = "metabo_file", 
        "Choose a metabolite target list", 
        placeholder = uploadMetaboFileName()
        )
    })
    
    output$demoMetaboTable <- renderDataTable({
     
      req(input$viewDemoMetabo, selTargetFilePath())

      if(input$viewDemoMetabo){
    
        target_path <- selTargetFilePath()
        if (stringr::str_ends(target_path, "csv")) {
          metabo_data<-metabolie.data.import.fun(target_path, IonPolarity = input$IonPolarity)
          datatable(metabo_data)
        } else if (stringr::str_ends(target_path, "library")) {
          metabo_data<-read_library_Fun(lib_dir = target_path, IonPolarity = input$IonPolarity)
         datatable(metabo_data)
        }
        
      }
      else{
        datatable(NULL)
      }

    })
    
    
    # show table check box
    
    output$metabo <- renderDataTable({
      
      req(input$showMetabo,selTargetFilePath())
      
      
      if (input$showMetabo && !is.null(selTargetFilePath())) {
        target_path<-selTargetFilePath()
        if (stringr::str_ends(target_path, "csv")) {
          metabo_data<-metabolie.data.import.fun(target_path, IonPolarity = input$IonPolarity)
        } else if (stringr::str_ends(target_path, "library")) {
          metabo_data<-read_library_Fun(lib_dir = target_path, IonPolarity = input$IonPolarity)
        }
        
        datatable( metabo_data )
      } else {
        
        datatable(NULL)
      }

    })
    

    
    # Download handlers for example csv file - PA14
    output$download5 <- downloadHandler(
      filename = function() {
        "Demo_PA14_target_list.csv"
      },
      content = function(file) {
        file.copy(demo_PA14Metabo_path, file)
      }
    )
    # Download handlers for example csv file - Si11
    output$download6 <- downloadHandler(
      filename = function() {
        "Demo_Si11_target_list.csv"
      },
      content = function(file) {
        file.copy(demo_Si11Metabo_path, file)
      }
    )


  observe({
    output$NLtable<-DT::renderDataTable({
      if(input$showNL){  
        get.NL()
         }
      else{
        NULL
        }
      })
    })
  observe({
    output$adductsTable<-DT::renderDataTable({
      if(input$showAdducts){  
        get.adduct()
      }
      else{
        NULL
      }
    })
  })
 

  output$NL_file <- renderUI({
    if (file.exists(default_NL_file_path)) {
      fileInput("NL_file", "Choose Neutral Loss File", placeholder = "Default file selected")
    } else {
      fileInput("NL_file", "Choose Neutral Loss File", placeholder = "No file selected")
    }
  })

  
  get.NL <- reactive({
  
    if (!is.null(input$NL_file) && input$NL_file$size > 0) {
        neutral_loss<-read.csv(input$NL_file$datapath, check.names = TRUE)
        }
      else if (file.exists(default_NL_file_path)) {
        neutral_loss<-read.csv(default_NL_file_path, check.names = TRUE)
      }
      if (!is.null(neutral_loss)) {
        if(any(input$IonPolarity%in%c("pos"))){
          neutral_loss_df<-neutral_loss[c("Accurate.Mass","Neutral.Loss","Pos")]
          
          NL_df<- neutral_loss_df[which(neutral_loss_df$Pos=="+"),]
          
          return(NL_df)
        }
        if(any(input$IonPolarity%in%c("neg"))){
          
          neutral_loss_df<-data.frame(neutral_loss[c("Accurate.Mass","Neutral.Loss","Neg")])
          
          NL_df<- neutral_loss_df[which(neutral_loss_df$Neg=="+"),]
          return(NL_df)
        }
      }
  })
  

  
  output$adduct_file <- renderUI({
   
      if (file.exists(default_adduct_file_path)) {
        fileInput("adduct_file","Choose adducts file", placeholder = "Default file selected")
        } else {
        fileInput("adduct_file","Choose adducts file", placeholder = "No file selected")
        }
    
  })
  
  get.adducts<- reactive({
    if (!is.null(input$adduct_file) && input$adduct_file$size > 0) {
        adduct_table<-read.csv(input$adduct_file$datapath, check.names = TRUE)
      }
    else if (file.exists(default_adduct_file_path)) {
      adduct_table<-read.csv(default_adduct_file_path, check.names = TRUE)
    }
    if ("Ion.mass" %in% colnames(adduct_table)) {
      # Check if any element in the "Ion.mass" column contains a space
      if (any(grepl(" ", adduct_table[["Ion.mass"]]))) {
        # Remove spaces from the "Ion.mass column
        adduct_table[["Ion.mass"]]<- gsub(" ", "", adduct_table[["Ion.mass"]])
         return(adduct_table)
        }
    }
    else{ cat("Column 'Ion.mass' not found in adduct_table\n")}
   
 }) 
    

  
  observe({
    output$adductsTable<-DT::renderDataTable({
      if(input$showAdducts){  
        get.adducts()
      }
      else{
        NULL
      }
    })
  })
  
  observe({
    if (input$IonPolarity == "pos") {
      updateSelectInput(session, "PIon",
                        choices = c("[M+H]+", "[M+Na]+"),
                        selected = c("[M+H]+", "[M+Na]+"))
    } else if (input$IonPolarity == "neg") {
      updateSelectInput(session, "PIon",
                        choices = c("[M-H]-"),
                        selected = "[M-H]-")
    }
  })
  
  get.PIon <- reactive({
    input$PIon
  })
  get.ppm <- reactive({
    input$ppm_exact
  })
  
  # match by PIon: [M+H]+
  get.match<-reactive({
    
    PIon<-get.PIon()
    ppm<-get.ppm()
    featureData<-get.df()
    
    if (is.null(selTargetFilePath())) {
      
      return(NULL)
      
    }
    
    compoundData<-get.metabo()
    
    # create RT NA column if RT column is missing
    if(!any("RT"==colnames(compoundData))){
      compoundData$RT <- rep(NA_real_, nrow(compoundData))
    }

    # import findMatch.cpp
    PI_res<-PImatch_fun(FT =featureData,Comp_data = compoundData,ppm = input$ppm_exact, 
                        PIon =input$PIon,diff_mz_thr =input$mz_diff_exact,diff_rt_thr=input$rt_thr_exact)

    PI_res_sub<-PI_res[which(sapply(PI_res, function(x) length(x$Feature_name)>0))]

    PI_match.df<- plyr::ldply(PI_res_sub, data.frame)
   
    return(PI_match.df)
  })
  
  # get adducts matching 
  find.adducts<-reactive({
    req(get.adducts())
    
    add_df<-get.adducts()
    
    if(input$IonPolarity=="pos"){
     
      add_df<-add_df[which(add_df$"Ion.polarity"=="pos"),]
      
      # adducts from MetaboCoreUtils library
      add_names<-adductNames(polarity = "positive")
      #keep what is unique in MetaboCoreUtils
      unique_add_names<-setdiff(add_names,add_df$Ion.name)
      
      unique_add_names_sub<-unique_add_names[which(!unique_add_names%in%c("[M+H+Na2]3+","[M+Na3]3+" ))]
    }
    if(input$IonPolarity=="neg"){
      # adduct table from the default path
      add_df<-add_df[which(add_df$"Ion.polarity"=="neg"),]
      
      # adducts from MetaboCoreUtils library
      add_names<-adductNames(polarity = "negative")
      #keep what is unique in MetaboCoreUtils
      unique_add_names<-setdiff(add_names,add_df$Ion.name)
      unique_add_names_sub<-unique_add_names[which(!unique_add_names%in%c("[M+H+Na2]3+","[M+Na3]3+" ))]
    }
    PIon<-get.PIon()
    ppm<-get.ppm()
    featureData<-get.df()
    
    if (is.null(selTargetFilePath())) {
     
      return(NULL)
      
    }
    compoundData<-get.metabo()
    
    n<-nrow(featureData)#
    # create rt NA column if it is missing
    if(!any("RT"==colnames(compoundData))){
       compoundData$RT <- rep(NA_real_, nrow(compoundData))
    }
   
    
    # import findMatch.cpp
    PI_res<-PImatch_fun(FT =featureData,Comp_data = compoundData,ppm = input$ppm_exact, PIon =input$PIon,diff_mz_thr =input$mz_diff_exact,diff_rt_thr=input$rt_thr_exact)
    # remove null lists
   
    PI_res_sub<-PI_res[which(sapply(PI_res, function(x) length(x$Feature_name)>0))]
    
    # generate adduct mass list
    adducts_mass_list<-list()
    for(i in seq_along(PI_res_sub)) {
      for(j in 1:length(PI_res_sub[[i]]$Feature_name)){
        M_val<-PI_res_sub[[i]][["M"]][j]
        adducts_mass_MetaboCoreUtils_temp<-mass2mz(M_val, adduct = unique_add_names_sub)
        adducts_mass<-c(adducts_mass_function(M = M_val,adduct_formula =add_df$Ion.mass),as.numeric(adducts_mass_MetaboCoreUtils_temp[1,]))
        adduct_name<-c(add_df$Ion.name,unique_add_names_sub)
        
        adducts_mass_list_temp <- list(Feature_name = PI_res_sub[[i]]$Feature_name[j],Feature_rt=PI_res_sub[[i]]$Feature_rt[j], 
                            Feature_mz = PI_res_sub[[i]]$Feature_mz[j],Comp_name = PI_res_sub[[i]]$Comp_name[j],
                            PI_name = PI_res_sub[[i]]$PI_name[j], M = PI_res_sub[[i]]$M[j],adducts_mass = adducts_mass, adduct_name = adduct_name)
        
        adducts_mass_list[[length(adducts_mass_list) + 1]] <- adducts_mass_list_temp
      }
    } 
      
   
    # search the feature table: find which feature is the adduct 
    adducts_anno_res<-adducts_anno_fun(PI_res = adducts_mass_list, FT = featureData,rt_thr = input$rt_thr_adducts,
                                        mz_thr = input$mz_diff_adducts,ppm = input$ppm_adducts)
    # mutate adduct and M as a new element 
    M_adducts_lists <- lapply(adducts_anno_res, function(x) {
      x$M_adducts <- paste(x$adduct_name,round(x$M, digits = 4))
      return(x)
    })
    # convert list to data frame and remove duplicated rows
    M_adducts_df<-unique(plyr::ldply(M_adducts_lists, data.frame))
    
    # keep adducts "feature_name" and " adducts_anno"
    M_adducts_df_sub<- M_adducts_df[c("feature_name","M_adducts")]
    
    PI_M_adducts<- M_adducts_df[c("PI_feature_name_vec","feature_name")]
   
    M_adducts_list<-list(M_adducts_df_sub,PI_M_adducts)
    
    
    return(M_adducts_list)
  })
  
  
  
  # Check C13 isotopes and multiply charge states
  
  get.isoCheck<-reactive({
    req(input$picker)
    req(selFilePath())
    # feature table order by mz 
    FT<-get.df()
    FT_sort<-FT[order(FT$mz),]

    # mutate a median intensity column for isotopes annotation
    FT_median_df<- FT_sort %>%
      dplyr::select(.,input$picker) %>% 
      dplyr::mutate_if(is.character,as.numeric) %>% 
      dplyr::mutate(median_value = apply(select(., all_of(input$picker)), 1, function(row) median(row, na.rm = TRUE))) %>% 
      data.frame()
    
    
    # Rcpp code- output 
    iso_anno<-Check_Iso_Charge(name =FT_sort$feature_name,mz = FT_sort$mz,
                               rt = FT_sort$rt,intensity =FT_median_df$median_value )
    
    FT_join<-data.frame(cbind(FT_sort[c("feature_name", "mz","rt")],iso_anno))
    
    return(FT_join)
  })
  
  iso_relation<-reactive({
    
    req(get.isoCheck())
    
    FT_join<-get.isoCheck()
   
    # split the [M]+ or [M+1]+ for identify the M and isotopes
    pattern <- "\\[(\\d+)\\]\\[(M\\+?)\\]?\\d*"
    
    FT_join_mutate<-FT_join %>% 
      mutate(iso_group=str_match(iso_anno, pattern)[,2],
             M_type=str_match(iso_anno, pattern)[,3]) %>% 
      na.omit() %>% 
      group_by(iso_group) %>% 
      data.frame()
    
    # Define a function to pair "M" with "M+" for a single iso_group to aasign the group index
    pair_m_mplus <- function(df) {
      m_df <- df %>% filter(M_type == "M") %>% select(feature_name)
      m_plus_df <- df %>% filter(M_type == "M+") %>% select(feature_name)
      
      paired_df <- expand.grid(m_df$feature_name, m_plus_df$feature_name, stringsAsFactors = FALSE)
      colnames(paired_df) <- c("PI", "ISF")
      paired_df_order<-paired_df[,c(2,1)]
      return(paired_df_order)
    }
    # Apply the pairing function to each iso_group
    output_df <- FT_join_mutate[c("feature_name","iso_group","M_type")] %>%
      group_by(iso_group) %>%
      group_split() %>%
      map_dfr(pair_m_mplus)
    
    return(output_df)
    
  })
  

  
  # find all pairs of features no correlation required
  find.all.NL<-reactive({
    req(get.adj_mat_NL())
    req(get.isoCheck())
    req(iso_relation())
    
    NL.df<-get.NL()
    adj_long.df<-get.adj_mat_NL()
  
    # import findNL_fun from cpp file
    find_NL<-findNL_fun(NL_data = NL.df,adj_long = adj_long.df,
                        diff_mz_thr = input$mz_diff_NL,
                        ppm=input$ppm_NL,diff_rt_thr = input$rt_thr_NL )
    find_NL.sub<-find_NL[which(sapply(find_NL, function(x) length(x$Feature_name1)>0))]
    find_NL.sub.df<-plyr::ldply(find_NL.sub, data.frame)
    
    return(find_NL.sub.df)
  })
  
  observe({
    
    output$full.NL= DT::renderDataTable(find.all.NL())
  })
  
  
  error_messages <- reactiveValues(messages = list())
 
  observe({
    req(input$picker)
    req(selFilePath())
    req(get.df())
    
    FT <- get.df()
    validate(
      need(ncol(df) > 0, "The uploaded file does not contain any columns.")
    )
    
    # check for non-numeric columns in the selected data
    non_numeric_cols <- sapply(FT[, input$picker, drop = FALSE], function(col) !is.numeric(col) )
    
    if (any(non_numeric_cols)) {

      non_numeric_col_names <- names(FT)[input$picker][non_numeric_cols]
      error_message <- paste("The following columns are not numeric and cannot be processed:", paste(non_numeric_col_names, collapse = ", "))
      error_messages$messages <- c(error_messages$messages, error_message)
      return()
    }
    
    tryCatch({
      mat <- FT %>%
        dplyr::select(., input$picker) %>%
        dplyr::mutate_if(is.character, as.numeric) %>% 
        mutate(across(where(is.numeric), ~na_if(., 0))) %>%  # Replace 0 with NA in numeric columns
        mutate(across(where(is.numeric), ~ifelse(is.na(.), NA_real_, log10(.)))) 
      
      if (!all(sapply(mat, is.numeric))) {
        stop("There are non-numeric columns.")
      }
    }, error = function(e) {
       e$message
    })
  })
  

  
  ## calculate adjacency matrix
  adj_matrix<-reactive({
    req(input$picker)
    req(get.df())
    
    FT<-get.df()
    
    mat<-FT %>% 
      dplyr::select(.,input$picker) %>%   # picked columns
      mutate(across(where(is.numeric), ~na_if(., 0))) %>%  # replace 0 with NA in numeric columns
      mutate(across(where(is.numeric), ~ifelse(is.na(.), NA_real_, log10(.)))) 
    

    t_mat<-t(mat)
    
    colnames(t_mat)<-get.df()$feature_name
    # point-to-point correlations, skipping NAs
    res.cor <- pairwiseCor(x = t_mat,method = input$cor_method) 

    
    colnames(res.cor)<-get.df()$feature_name
    rownames(res.cor)<-get.df()$feature_name
    # get rid of low correlations and diagnal values
    res.cor[res.cor < input$cor_thr ] <- NA
    # Set the lower triangle, including the diagonal, to NA
    res.cor[lower.tri(res.cor, diag = FALSE)] <- NA
    return(res.cor)
    
  })
  
  
  get.adj_mat<-reactive({
    req(get.df())
    req(adj_matrix())
    adj_mat.long<-na.omit(data.frame(as.table(adj_matrix())))
    
    # left join to add Var1, Var2,Freq,mz_x,mz_y,rt_x,rt_y,rt_diff,mz_diff
    adj.full<-left_join_and_mutate_fun(FT_df=get.df(),long_df = adj_mat.long,rt_thr=input$rt_thr)
  

    return(adj.full)
  })
  get.adj_mat_NL<-reactive({
    
    req(get.df())
    req(get.isoCheck())
    req(iso_relation())
   
    iso_FT_name<-iso_relation()$ISF
    FT<-get.df()
    FT_deiso<-FT[which(!FT$feature_name%in%iso_FT_name),]
    
    combs <- combn(FT_deiso$feature_name, 2)
    
    combinations_df <- as.data.frame(t(combs), stringsAsFactors = FALSE)
    names(combinations_df) <- c("Var1", "Var2")
    # left join to add Var1, Var2,mz_x,mz_y,rt_x,rt_y,rt_diff,mz_diff, annotation using Rcpp code
    adj.full<-left_join_and_mutate_NL_fun(FT_df=FT_deiso,long_df = combinations_df,rt_thr=input$rt_thr)
  
    
    return(adj.full)
  })
  #####################################################################################
  # MS2 and PI annotation
  # input for assign group index
  # output: data frame with two columns   
  
  get.ms2_PI<-reactive({

    file_path <- selectedMS2FilePath()
    if (is.null(file_path) || file_path == "" || !file.exists(file_path)) {

      
      # Return an empty data frame with columns "ISF" and "PI"
      return(data.frame(ISF = " ", PI = " "))
    }
    
    
    if(any(grepl("mzXML",selectedMS2FilePath())|grepl("mzML",selectedMS2FilePath()))){
      req(selectedMS2FilePath())
      req(get.ms2_precursor())
      req(get.df())
      # feature table order by mz 
      FT<-get.df()
     
      mzMatch_RES_df<-get.ms2_precursor()[[1]]

      
      # input for cpp:assign the group index
      MS2_PI_input<-data.frame(mzMatch_RES_df[c("feature_name","Precursor_matched_FT_name")])
      colnames(MS2_PI_input)<-c("ISF","PI")

      return(MS2_PI_input)
      
    }
    
    if(any(grepl("mgf",selectedMS2FilePath()))){
      
      req(get.ms2_precursor_MZmine())
      req(get.df())
      # feature table order by mz 
      FT<-get.df()
      FT_sort<-FT[order(FT$mz),]
      mzMatch_RES_df<-get.ms2_precursor_MZmine()[[1]]
      

      # input for cpp:assign the group index
      MS2_PI_input<-data.frame(mzMatch_RES_df[c("ms2_feature_name","PI_feature_name")])
      colnames(MS2_PI_input)<-c("ISF","PI")

      return(MS2_PI_input)
      
    }
    
  })
  
  # NL data frame as input
  get.NL.input<-reactive({
    req(find.all.NL())
    
    # input for cpp:assign the group index
    NL_PI_input<-data.frame(find.all.NL()[c("Feature_name1","Feature_name2")])
    colnames(NL_PI_input)<-c("PI","ISF")
    NL_PI_input_reorder<-NL_PI_input[,c(2,1)]

    return(NL_PI_input_reorder)
    
  })
  
  
  # join all NL, MS2, adduct and PI data frame 
  
  join_NL_MS2_PI<-reactive({
    req(get.ms2_PI())
    req(get.NL.input())
    #req(find.adducts())
    req(iso_relation())
    
    if (!is.null(selTargetFilePath())) {
      # change column name
      PI_adducts_df<-find.adducts()[[2]]
      colnames(PI_adducts_df)<-c("PI","ISF")
      # relocate the order of PI_adducts_df
      PI_adducts_df_order<-PI_adducts_df[, c(2, 1)]
    
      join_df<-unique(rbind(get.ms2_PI(),get.NL.input(),PI_adducts_df_order,iso_relation()))
   
      output_group<- groupFeatures(join_df)
   
      return(output_group)
    }
    else{

      join_df<-unique(rbind(get.ms2_PI(),get.NL.input(),iso_relation()))
      
      output_group<- groupFeatures(join_df)
      
      return(output_group)
    }
  })
  
  #######################################################################################################
  # assign correlation group
  get_corGroup<-reactive({
    req(get.adj_mat())
    
    get.adj_matrix<-get.adj_mat()
    # order by Freq
    get.adj_order<-get.adj_matrix[order(get.adj_matrix$Freq,decreasing = T),]
    # keep the unique
   
    get.adj_order_unique<-get.adj_order[which(get.adj_order$Var1!=get.adj_order$Var2),]
    # assign group index
    groupCorFeatures_res<-groupCorFeatures(adj_long =get.adj_order_unique,threshold=input$cor_thr)
    return(groupCorFeatures_res)
  })
  
  #####################################################################################
  ## match feature table mz to precursor mz from MS2 file
  
  get.ms2_precursor<-reactive({
    req(selectedMS2FilePath())
    req(selFilePath())
    req(get.isoCheck())
    req(iso_relation())   

    FT<-get.df()

    # check if the file exists
    if (is.null(selectedMS2FilePath())) {
      message("Selected MS2 file path does not exist. Skipping MS2 processing.")
      return(NULL) 
     
    }
    ms2 <-MSnbase::readMSData(files = selectedMS2FilePath(), msLevel. = 2, mode = "onDisk")
    

    # match the precursor mz to feature table mz-cpp code
    MS2match_res<-MS2match(MS_OBJ = ms2,FT =FT, 
                             mz_diff_precursor = input$mz_diff_precursor,
                             rt_thr_precursor = input$rt_thr_precursor,ppm = input$ppm_precursor)
    
    MS2match_res_sub<-MS2match_res[which(sapply(MS2match_res, function(x) length(x$Feature_name)>0))]

    # keep the Precursor which has the maximum intensity
    get_maxi_list <- lapply(MS2match_res_sub, function(x) {
        index <- which.max(x$Precursor_intensity)
        output <- list(
          Feature_name = x$Feature_name[index],
          Feature_mz = x$Feature_mz[index],
          Feature_RT = x$Feature_RT[index],
          Precursor_mz = x$Precursor_mz[index],
          Precursor_rt = x$Precursor_rt[index],
          Precursor_intensity = x$Precursor_intensity[index],
          mz_Vals = x$mz_Vals[index]
        )
        return(output)
      })
      
      mzMatch_RES<-MS2_mzMatch(input_list = get_maxi_list,FT = FT,
                               mz_diff_MS2 = input$mz_diff_MS2,
                               rt_thr_MS2 = input$rt_thr_MS2,ppm = input$ppm_MS2)
      
      mzMatch_RES_df<-plyr::ldply(mzMatch_RES, data.frame)

      ISF_anno<-ISFAnno_fun(mzMatch_df =mzMatch_RES_df,FT = FT)
     
      output_df<-data.frame(cbind(FT[c("feature_name","mz","rt")],ISF_anno))

      output_list<-list(mzMatch_RES_df,output_df)
      return(output_list)
  
  })
  # check if the correlation of intensity condition holds: in-source fragment match 
  
  get.ms2_precursor_MZmine<-reactive({
    
    req(selectedMS2FilePath())
    req(selFilePath())
    
    # feature table order by mz 
    FT<-get.df()
    FT_sort<-FT[order(FT$mz),]
    
    if (is.null(selectedMS2FilePath())) {
      message("Selected MS2 file path does not exist. Skipping MS2 processing.")
      return(NULL)  # Skip the process if the file does not exist
      
    }
      
    MS2_mgf_import<-MS2_mgf_import_fun(mgf_file=selectedMS2FilePath())
      
    ISFMZmine_list<-ISFMZmine_fun(MZmine_list =MS2_mgf_import,FT = FT_sort,rt_thr_MS2 = input$rt_thr_MS2,
                                    ppm = input$ppm_MS2,ms2_mz_diff = input$mz_diff_MS2 )
      
      
    ISFMZmine_df<-plyr::ldply(ISFMZmine_list, data.frame)
    # assign: PI match and MS2 match to feature table
    ISF_anno<-ISFMZmine_assign_fun(PI_MS2_df=ISFMZmine_df,FT=FT_sort)  #mzMatch_df_join_max
    
    # cbind FT name, mz, rt, ISF annotation
    ISF_anno_mzmine_df<-data.frame(cbind(FT_sort[c("feature_name","mz","rt")],ISF_anno))
      
    #progress$inc(1/1, detail = paste("Doing part", 1))
 
    ISF_anno_mzmine_list<-list(ISFMZmine_df,ISF_anno_mzmine_df) 
    return(ISF_anno_mzmine_list)
 
    
  })

  ######################################################################################
  # subset of NL - select row from the output feature table
  find.selected.NL<-reactive({
    req(input$Output_FT_rows_selected)
    req(get.FT.output())
    req(get.df())
    req(find.all.NL())
    req(selectedFeatureName())
    
    sel_fea.name<-selectedFeatureName()
    NL.df<-get.NL()
    adj_long.df<-find.all.NL()
    fea_table<-get.df()
    compoundData<-get.metabo()
    PIon<-get.PIon()

    # selected feature names
    adj_long.df.sel<-adj_long.df[which(adj_long.df$Feature_name1%in%sel_fea.name|adj_long.df$Feature_name2%in%sel_fea.name),]

    return(adj_long.df.sel)
  })
 
  observe({
    req(input$FT_rows_selected)
    output$selected.NL= DT::renderDataTable(find.selected.NL())
  })
  

 
  
  output$picker <- renderUI({
    req(selFilePath())
    filter_out_column_names <- c("feature_name", "compound", "mz", "mzmin", "mzmax", "rt")
    df <- get.df()
    validate(
      need(ncol(df) > 0, "The uploaded file does not contain any columns.")
    )
   
    pattern_columns <- c()
    
    if (any(grepl("_[A-Z]$", colnames(df)))) {
     
      pattern_columns <- colnames(df)[grepl("_[A-Z]$", colnames(df))]
      
    } else if (any(grepl(".mzXML.Peak.area", colnames(df)))) {
     
      pattern_columns <- colnames(df)[grepl(".mzXML.Peak.area", colnames(df))]
      
    } else if (any(grepl("_(\\d+)$", colnames(df)))) {
     
      pattern_columns <- colnames(df)[grepl("_(\\d+)$", colnames(df))]
    } else {
      
      pattern_columns <- colnames(df)[which(!colnames(df)%in%c("featureidx","CV", "pvalue","qvalue","feature_name","npeaks", "compound", "mz", "mzmin", "mzmax", "rt","rtmin" , "rtmax" , "npeaks" ))]

    }
   
    
    if (any(grepl("Peak.area", colnames(df)))) {
   
      pickerInput(inputId = 'picker',
                  label = 'Select columns ',
                  choices = colnames(df),
                  selected = colnames(df)[which(grepl("Peak.area", colnames(df)))],
                  multiple = TRUE)
    } else {
    
      pickerInput(inputId = 'picker',
                  label = 'Select columns ',
                  choices = colnames(df),
                  selected = pattern_columns,
                  multiple = TRUE)  
    }
  })
  
  observe({
    output$find.metabo <- DT::renderDataTable(get.match())
  })
  
  observe({
    output$Iso_check_table <- DT::renderDataTable(get.isoCheck())
  })
  
  ### VisNetwork network plot
  
  
  get.ft<-reactive({
    FT<-get.df()
    return(FT)
  })
  
  
  ############################################################################################
  # the feature table for output
  get.FT.output <- reactiveVal(NULL)
  
  observeEvent(input$run, {
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive
      on.exit(progress$close())
      progress$set(message = "Start calculation...")      

     
      Sys.sleep(2) 
      
      progress$set(message = "Feature annotation", value = 0)
      
   
      # feature table order by mz 
      FT<-get.df()
      fea.table<-FT[order(FT$mz),]
     
      # check if the target list is uploaded
      if (is.null(selTargetFilePath()) ){
        fea_anno.output <- data.frame(feature_name = FT[,"feature_name"]) %>% 
          mutate(metabolite_annotation = rep("", nrow(FT))) %>% 
          data.frame()
        
        progress$inc(1/7, detail = paste("No target list file- Skipping metabolite annotation step", 1))
      }
      
      else{
        MH.df<-get.match()# mutate a column by paste compound names and PIon names
        progress$inc(1/7, detail = paste("Metabolite identification", 1))
        Sys.sleep(2) #
      
        MH.df_sub <- MH.df%>%
          dplyr::mutate(Comp_PI_name = paste(Comp_name, PI_name, sep = " ")) %>% 
          dplyr::select(Feature_name,Comp_PI_name)
      
      
        # self defined function in helper.r: convert long to wide format
        MH.df.output<-long_to_wide.fun(df=MH.df_sub,col_name1 =Feature_name,col_name2=Comp_PI_name )
        rm(MH.df)
     
      
        fea_anno.df<-MH.df.output

        fea_anno.output<-long_to_wide.fun(df=fea_anno.df,col_name1 =feature_name,col_name2=feature_annotation )
        colnames(fea_anno.output)[2]<-"metabolite_annotation"
      }

      # check if any mzXML, mzML, or mgf file  in selectedMS2FilePath()
      if (!is.null(selectedMS2FilePath())) {
        if (any(grepl("mzXML", selectedMS2FilePath()) | grepl("mzML", selectedMS2FilePath()))) {
          req(get.ms2_precursor())
          
          output_ms2_precursor_table <- get.ms2_precursor()[[2]]
          progress$inc(1/7, detail = paste("MS2 match annotation", 2))
          Sys.sleep(2)
          
        }
        
        if (any(grepl("mgf", selectedMS2FilePath()))) {
          req(get.ms2_precursor_MZmine())
          output_ms2_precursor_table <- get.ms2_precursor_MZmine()[[2]]
        }
        
        progress$inc(1/7, detail = paste("MS2 match annotation", 2))
        Sys.sleep(2)
        
      } else {
      
        output_ms2_precursor_table <- data.frame(feature_name = FT[,"feature_name"]) %>% 
          mutate(ISF_anno = rep("", nrow(FT))) %>% 
          data.frame()
        
        progress$inc(1/7, detail = paste("No MS2 files found. Skipping MS2 annotation step", 2))
        Sys.sleep(2)
      }
      
      

      output_Iso_check_table<-get.isoCheck() 
      
      progress$inc(1/7, detail = paste("Isotopes annotation", 3))
      Sys.sleep(2) #
      
      # adducts annotation data frame
      if (is.null(selTargetFilePath())) {
        M_adducts_df_merge <- data.frame(feature_name = FT[,"feature_name"]) %>% 
          mutate(adducts_anno = rep("", nrow(FT))) %>% 
          data.frame()
        progress$inc(1/7, detail = paste("Skipping adducts annotation", 4))
        Sys.sleep(2) #    
       
      }
      else{      
      long_output_adducts<-unique(find.adducts()[[1]])
      # long to wide: adducts output
      M_adducts_df_merge<-long_to_wide.fun(df = long_output_adducts,col_name1 =feature_name, col_name2 = M_adducts)
      # rename the column 
      colnames(M_adducts_df_merge)[2] <- "adducts_anno"  
      progress$inc(1/7, detail = paste("Adducts annotation", 4))
      Sys.sleep(2) #
      }


      NL.df<-find.all.NL()

      NL.df$fea_anno<-apply(NL.df, 1,function(x) paste0(paste0(x["Feature_name1"],"-"),x["NL_Names"]))
      
      NL.fea2.output<-long_to_wide.fun(df=NL.df,col_name1 =Feature_name2,col_name2=fea_anno )
      NL_anno.output<-long_to_wide.fun(df=NL.fea2.output,col_name1 =feature_name,col_name2=feature_annotation )
      colnames(NL_anno.output)[2]<-"neutral_loss_annotation"
      
      rm(NL.df)
      progress$inc(1/7, detail = paste("Neutral loss annotation", 5))
      Sys.sleep(2) 
      
      # assign group index
      group_index<-join_NL_MS2_PI()
      
      progress$inc(1/7, detail = paste("assign group index", 6))
      Sys.sleep(2) #
      # correlation groups
      group_Cor_index<-get_corGroup()[[1]]%>%
        group_by(cor_group) %>%
        filter(n() > 1) %>%
        ungroup() %>% 
        data.frame()
      
      progress$inc(1/7, detail = paste("assign correlation group index", 7))
      Sys.sleep(2) #
      # multiple data frame mergeing
      list_of_dfs <- list(fea.table[c("feature_name","mz","rt")], fea_anno.output,NL_anno.output,output_ms2_precursor_table[c("feature_name","ISF_anno")],M_adducts_df_merge,
                          output_Iso_check_table[c("feature_name","iso_anno")], group_index,group_Cor_index)
      
      # Use reduce to left_join all data frames by "feature_name"
      merged_df <- purrr::reduce(list_of_dfs, left_join, by = "feature_name")
    
      rm(M_adducts_df_merge)
      rm(output_Iso_check_table)
      rm(group_index)
      rm(group_Cor_index)
      
      # output table format
      output<-merged_df %>% 
        dplyr::select("feature_name","mz","rt","metabolite_annotation","adducts_anno",
                      "neutral_loss_annotation","ISF_anno","iso_anno",
                      "group","cor_group") %>% 
        dplyr::mutate(mz = round(mz, 4),
                      rt=round(rt, 1),
                      user_anno = rep("",nrow(fea.table))) %>% 
        data.frame()
      
      get.FT.output(output)
   
  })

  

  selectedFeatureName <- reactiveVal()
  
  observe({
    df <- get.FT.output()
    groupNumbers <- sort(as.numeric(str_extract(unique(df$group), "\\d+")))
    corgroupNumbers <- sort(as.numeric(str_extract(unique(df$cor_group), "\\d+")))
    
    groups <- paste0("group", groupNumbers)
    corgroups <- paste0("corgroup", corgroupNumbers)
    
    choices <- c("All" = "All", groups, corgroups)
    updateSelectInput(session, "combinedFilter", choices = choices)
   
  })
  
  filteredData <- reactive({
    df <- get.FT.output()  # Assuming this is your base data fetching function
    if (input$combinedFilter == "All") {
      df
    } else if (input$combinedFilter %in% df$group) {
      df %>% filter(group == input$combinedFilter)
    } else if (input$combinedFilter %in% df$cor_group) {
      df %>% filter(cor_group == input$combinedFilter)
    } else {
      df
    }
  })
  output$Output_FT <- renderDT({
    datatable(filteredData(), editable = TRUE, selection = 'single', options = list(
      pageLength = 10,
      autoWidth = TRUE,
      searchHighlight = TRUE,
      searching = TRUE
    ), filter = 'top')
  })

  observeEvent(input$Output_FT_rows_selected, {
    selectedRow <- input$Output_FT_rows_selected
    req(selectedRow)  # Ensure a row is selected
    df_current <- filteredData()  # Use the filtered dataset
    if (!is.null(selectedRow) && selectedRow > 0 && selectedRow <= nrow(df_current)) {
      selectedFeatureName(df_current[selectedRow, "feature_name", drop = TRUE])

    }
  })

  # box plot with selected columns and selected rows in the output feature table
  box_df<-reactive({
    req(input$Output_FT_rows_selected)
    req(input$picker)
    req(input$Output_FT_rows_selected)
    req(get.plot.output())
    req(selectedFeatureName())
    if (is.null(input$picker) || length(input$picker) == 0) return(NULL)
    
    if (any(grepl("_[A-Z]$", input$picker))) {
      # any column ends with _A, _B, _C, ... _Z
      col_names <- sub("(.*)_[A-Z]$", "\\1", input$picker)
      
    } else if (any(grepl(".mzXML.Peak.area", input$picker))) {
      col_names <- gsub(".mzXML.Peak.area", "", input$picker)
      
    } else if (any(grepl("_(\\d+)$", input$picker))) {
      # any column ends with _1, _2, _3, ... _n
      col_names <- sub("(.*)_\\d+$", "\\1", input$picker)
      
    } else {
      col_names <- input$picker
    }
    
    FT<-get.df()
  
    g.df<-get.plot.output()
    Visg_ft<-toVisNetworkData(g.df)
    
    nodes_df<-Visg_ft$nodes
    nodes_sub <-nodes_df
    
    FT_name<-nodes_sub$id 

    sel_df<-FT[which(FT$feature_name%in%FT_name),c("feature_name",input$picker)]
    
    t_df<-setNames(data.frame(t(sel_df[,-1])), sel_df[,1])
    df_wide <- tibble::rownames_to_column(t_df, "Sample")
    df_wide$sample_group<-as.factor(col_names)
    
    df_long<-tidyr::pivot_longer(df_wide,cols = starts_with("FT"), names_to = "FT_name", values_to = "y_value") %>% data.frame()
    
    df_long$FT_name<-as.factor(df_long$FT_name)
   
    df_long[df_long == 0] <- NA
    
    return(df_long)
  })
  



 
  plotHeight <- reactive({
    if(input$boxPlotCardMaximized) {
      "800px"
    } else {
      "500px"
    }
  })
  

  output$boxplot <- renderUI({
    plotOutput("box_Plot", height = plotHeight(), width = "100%")
  })

  output$box_Plot <- renderPlot({
    req(box_df())
    df_long <- box_df()
    df_filtered <- df_long %>%
      filter(!is.na(y_value) & y_value > 0)

    # changes in maximization state or window resize
    input$boxPlotCardMaximized
    plot_width <- session$clientData$output_boxplot_width


    gg <- ggplot(df_filtered, aes(x=sample_group, y=log10(y_value), fill=sample_group)) +
      geom_boxplot() +
      ylab("Log10 intensity") +
      facet_wrap(~FT_name) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    vals$gg <- gg
    print(gg)
  }, width = "auto")

  
  network_plotHeight <- reactive({
    if(input$networkPlotCardMaximized) {
      "800px"
    } else {
      "500px"
    }
  })
  
  output$output_network <- renderUI({
    visNetworkOutput("output_network_plot", height = network_plotHeight(), width = "100%")
  })
  # define a reactive value to store the visNetwork plot
  networkPlotReactive <- reactive({
    req(get.nodes())
    req(get.edges())
    
    nodes_df<-get.nodes()
    edges_df<-get.edges()
    visNetwork(nodes = nodes_df, edges = edges_df) %>% 
      visIgraphLayout(layout = "layout_in_circle") 
  })
  # create the visNetwork plot for output
  output$output_network_plot <- renderVisNetwork({
    networkPlotReactive()  # Use the reactive expression here
  })
  

  get.plot.output<-reactive({

    req(get.df())
    req(get.FT.output())
    req(selectedFeatureName())
    req(get_corGroup())
    
    
    
    sel_fea.name<- selectedFeatureName()

    adj_mat.list<-get_corGroup()[[2]]
  
  
    find_in_df <- function(df, search_FT) {
      any(df$Var1 == search_FT | df$Var2 == search_FT)
    }
    
    search_value <- sel_fea.name
    
    FT_index <- which(sapply(adj_mat.list, find_in_df, search_FT = search_value)==TRUE)
    if(length(FT_index)>0)
      {
      adj_mat.long<-adj_mat.list[[FT_index]]
      adj.full<-left_join_and_mutate_fun(FT_df=get.df(),long_df = adj_mat.long,rt_thr=input$rt_thr)
      g_ft<- graph_from_data_frame(adj.full, directed=FALSE)
      return(g_ft)
    }
    else{
      adj_mat.long<-data.frame(Var1=search_value,Var2=search_value,Freq=1)
      adj.full<-left_join_and_mutate_fun(FT_df=get.df(),long_df = adj_mat.long,rt_thr=input$rt_thr)
     
      g_ft<- graph_from_data_frame(adj.full, directed=FALSE)
     
      return(g_ft)
    }
      

  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("Feature_table_output", Sys.Date(), ".csv")
    },
    content = function(file) {
      accumulated_data <- modified_data() 
    
      if (is.null(accumulated_data)) {
        accumulated_data <- get.FT.output() # 
      }
      tryCatch({
        write.csv(accumulated_data, file, row.names = FALSE)
      }, error = function(e) {
        print(paste("Error writing CSV:", e))
      })
    }
  )
  vals <- reactiveValues()

  output$downloadboxPlot <- downloadHandler(
    filename = function() { paste0("Box_Plot_", Sys.Date(), ".png") },  
    content = function(file) {
      ggsave(file, plot = vals$gg, device = "png", width = 10, height = 8, units = "in", dpi = 300)
    }
  )
  
  get.nodes<-reactive({
    
    FT<-get.df()
    
    g.df<-get.plot.output()
    
    NL.df<-find.selected.NL()
    if (!is.null(selTargetFilePath())){
      MH.df<-get.match()
      
      isoCheck.df <- get.isoCheck()
      
      
      Visg_ft<-toVisNetworkData(g.df)
      nodes_df<-Visg_ft$nodes
      nodes_sub <-nodes_df
      
      find_mz<-paste0("m/z: ",round(FT$mz[match(nodes_df$id,FT$feature_name)],4))
      
      find_rt<-paste0("rt: ",round(FT$rt[match(nodes_df$id,FT$feature_name)],1))
      
      # label nodes are [M+H]+ PI found and isotopes
      
      find_mh<-MH.df$Feature_name
      
      # label nodes are adducts found
      find_adduct_df<-unique(find.adducts()[[1]])
      
      find_adduct_name<-find_adduct_df$feature_name
      
      # add text on node
      nodes_sub$color<-ifelse(nodes_sub$id%in%c(find_mh,find_adduct_name),"#86b9ff","#D2E5FF")
      
      # add text on node title
      nodes_sub$find.metabo<-rep("",nrow(nodes_sub))

      concat_comp_names <- function(id) {
        comp_names <- MH.df$Comp_name[MH.df$Feature_name == id]
        paste(comp_names, collapse = '<br/>')
      }
      # concatenate add_names for each id
      concat_add_names <- function(id) {
        add_names <- find_adduct_df$M_adducts[find_adduct_df$feature_name == id]
        paste(add_names, collapse = '<br/>')
      }
      nodes_sub$find.metabo <- purrr::map_chr(nodes_sub$id, concat_comp_names)
      nodes_sub$find.add <- purrr::map_chr(nodes_sub$id, concat_add_names)
      
      iso_feat_name<-isoCheck.df[grep("\\[M\\+",isoCheck.df$iso_anno),"feature_name"]
     
      iso<-FT$feature_name[which(FT$feature_name%in%iso_feat_name)]
      nodes_sub$color[which(nodes_sub$id%in%iso)]<-"#D3D3D3"
      
      ### node title for HTML
      nodes_sub$shape<-rep("box",nrow(nodes_sub))
      nodes_sub$sub_title<-paste0(find_mz,paste0("<br>",find_rt))
      nodes_sub$label<-find_mz
      nodes_sub$title<-paste0(nodes_df$id,"<br>",nodes_sub$sub_title,"<br>",nodes_sub$find.add,"<br>",nodes_sub$find.metabo)
      return(nodes_sub)
    }
    
    else{

      isoCheck.df <- get.isoCheck()
      
      Visg_ft<-toVisNetworkData(g.df)
      nodes_df<-Visg_ft$nodes
      nodes_sub <-nodes_df
      
      find_mz<-paste0("m/z: ",round(FT$mz[match(nodes_df$id,FT$feature_name)],4))
      
      find_rt<-paste0("rt: ",round(FT$rt[match(nodes_df$id,FT$feature_name)],1))

      find_adduct_df<-unique(find.adducts()[[1]])
      
      find_adduct_name<-find_adduct_df$feature_name
      
      # add text on node
      nodes_sub$color<-rep("#D2E5FF",nrow(nodes_sub)) 
      
      # add text on node title
      
      nodes_sub$find.metabo<-rep("",nrow(nodes_sub))
      
      iso_feat_name<-isoCheck.df[grep("\\[M\\+",isoCheck.df$iso_anno),"feature_name"]
  
      iso<-FT$feature_name[which(FT$feature_name%in%iso_feat_name)]
      nodes_sub$color[which(nodes_sub$id%in%iso)]<-"#D3D3D3"
      
      ### node title for HTML
      nodes_sub$shape<-rep("box",nrow(nodes_sub))
      nodes_sub$sub_title<-paste0(find_mz,paste0("<br>",find_rt))
      nodes_sub$label<-find_mz
      nodes_sub$title<-paste0(nodes_df$id,"<br>",nodes_sub$sub_title,"<br>",nodes_sub$find.add,"<br>",nodes_sub$find.metabo)
      return(nodes_sub)
    }
    
    
    })

  get.edges <- reactive({

    FT <- get.df()
    g.df <- get.plot.output()
    NL.df <- find.selected.NL()
    
    if (!is.null(selTargetFilePath())){
      
      MH.df <- get.match()
      isoCheck.df <- get.isoCheck()
      
      Visg_ft <- toVisNetworkData(g.df)
      nodes_df <- Visg_ft$nodes
      edges <- Visg_ft$edges

      edges_sub <- edges[edges$from != edges$to, ]
      
      if (nrow(edges_sub) > 0) {
        
        cor_r <- paste0("r: ", round(edges_sub$Freq, 2))
        delta_mz <- paste0("&Delta;", "m/z: ", round(abs(edges_sub$mz_x - edges_sub$mz_y), 2))
        
        concat_NL_names <- function(from, to) {
          matches <- NL.df %>%
            filter((Feature_name1 == from & Feature_name2 == to) | (Feature_name1 == to & Feature_name2 == from))
          paste(matches$NL_Names, collapse = '<br/>')
        }
        
        edges_sub$find.NL.text <- purrr::pmap_chr(list(edges_sub$from, edges_sub$to), concat_NL_names)
        
        edges_sub$title <- ifelse(edges_sub$find.NL.text != "", 
                                  paste0(cor_r, "<br>", delta_mz, "<br>", edges_sub$find.NL.text),
                                  paste0(cor_r, "<br>", delta_mz, "<br>"))
        
        if (nrow(edges_sub) > 1) {
          edges_sub$width <- rep(1, nrow(edges_sub))
          edges_sub$width[which(edges_sub$find.NL.text != "")] <- 3
        }
        
      }
      
      return(edges_sub)
      
    }
    else{
      #MH.df <- get.match()
      isoCheck.df <- get.isoCheck()
      
      Visg_ft <- toVisNetworkData(g.df)
      nodes_df <- Visg_ft$nodes
      edges <- Visg_ft$edges
      
      edges_sub <- edges[edges$from != edges$to, ]
      
      if (nrow(edges_sub) > 0) {
        
        cor_r <- paste0("r: ", round(edges_sub$Freq, 2))
        delta_mz <- paste0("&Delta;", "m/z: ", round(abs(edges_sub$mz_x - edges_sub$mz_y), 2))
        
        concat_NL_names <- function(from, to) {
          matches <- NL.df %>%
            filter((Feature_name1 == from & Feature_name2 == to) | (Feature_name1 == to & Feature_name2 == from))
          paste(matches$NL_Names, collapse = '<br/>')
        }
        
        edges_sub$find.NL.text <- purrr::pmap_chr(list(edges_sub$from, edges_sub$to), concat_NL_names)
        
        edges_sub$title <- ifelse(edges_sub$find.NL.text != "", 
                                  paste0(cor_r, "<br>", delta_mz, "<br>", edges_sub$find.NL.text),
                                  paste0(cor_r, "<br>", delta_mz, "<br>"))
        
        if (nrow(edges_sub) > 1) {
          edges_sub$width <- rep(1, nrow(edges_sub))
          edges_sub$width[which(edges_sub$find.NL.text != "")] <- 3
        }
        
      }

      return(edges_sub)
    }
    
  })
  
  
  
  # Update selected feature name based on row selection
  observeEvent(input$Output_FT_rows_selected, {
    selectedRow <- input$Output_FT_rows_selected
    req(selectedRow)  
    df_current <- filteredData()  
    if (!is.null(selectedRow) && selectedRow > 0 && selectedRow <= nrow(df_current)) {
      selectedFeatureName(df_current[selectedRow, "feature_name", drop = TRUE])

    }
  })



  output$downloadNetworkPlot <- downloadHandler(
    filename = function() {
      paste0("network_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot <- networkPlotReactive()
      
      if (inherits(plot, "htmlwidget")) {
        tempHtml <- tempfile(fileext = ".html")
        saveWidget(plot, tempHtml, selfcontained = TRUE)
        
        b <- ChromoteSession$new()
        tempHtmlUrl <- paste0("file://", normalizePath(tempHtml))
        
        b$Page$navigate(tempHtmlUrl)
        Sys.sleep(2) 
        
        b$Emulation$setDeviceMetricsOverride(width = 1000, height = 600, deviceScaleFactor = 0, mobile = FALSE)
        screenshot <- b$Page$captureScreenshot()
        b$close()
        
        tempImage <- tempfile(fileext = ".png")
        writeBin(base64enc::base64decode(screenshot$data), tempImage)
        
       
        image <- image_read(tempImage)
        image_cropped <- image_crop(image, "600x600+200+0")
        image_write(image_cropped, path = file) 
        
        unlink(tempHtml)  
        unlink(tempImage)
      } else {
        stop("Plot is not a visNetwork (htmlwidget) object.")
      }
    }
  )
  
  modified_data <- reactiveVal(NULL)
  
  observeEvent(input$Output_FT_cell_edit, {
    req(get.FT.output())
    
    info <- input$Output_FT_cell_edit
    
    if (!is.null(get.FT.output()) && !is.null(info$value)) {
      
      current_modified_data <- modified_data()
      
      if (is.null(current_modified_data)) {
        current_modified_data <- get.FT.output()
      }
      
      current_modified_data[info$row, info$col] <- info$value
      
      modified_data(current_modified_data)
    }
  })

}


