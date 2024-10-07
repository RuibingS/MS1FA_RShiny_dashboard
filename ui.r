library(shiny)
library(bs4Dash)
library(visNetwork)
library(DT)
library(shinythemes) 
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)


input.cor_methods<-c("pearson", "kendall", "spearman")


ui <- dashboardPage(

  header = dashboardHeader(
    title = dashboardBrand(
      title = tagList(
        span("MS1FA", style = "font-weight: bold;")
      ),
    color = "primary"),
    skin = "light",
    status = "white",
    border = TRUE,
    sidebarIcon = icon("bars"),
    controlbarIcon = icon("th"),
    fixed = FALSE
   ),
  
  sidebar = dashboardSidebar(
    skin = "light",
    status = "primary",
    elevation = 3,
    # Custom CSS to change tab background and text color
    tags$style(HTML('
      .sidebar .nav-sidebar .nav-item .nav-link {
        background-color: white !important;
        color: black !important;
      }
    ')),
    class="sidebar",
    sidebarMenu(
      sidebarHeader(""),
      menuItem(
        "Files upload",
        tabName = "files",
        icon = icon("file-upload")
      ),
      menuItem(
        "Select parameters",
        tabName = "parameters",
        icon =icon("sliders") 
      ),
      menuItem(
        "Output",
        tabName = "output",
        icon =icon("table") 
      )
    )
  ),
  body = dashboardBody(    
    useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "www/styles.css")),
    tags$head(
      tags$style(HTML("
        @media print {
            .main-sidebar { 
                display: block !important;
                position: fixed;
                top: 0;
                left: 0;
                width: 250px; /* Adjust based on your sidebar width */
                height: 100%;
                z-index: 1000;
            }
            .content-wrapper {
                margin-left: 250px; /* Should match the sidebar width */
                padding: 20px; /* Optional: Add padding for better readability */
            }
            .main-header, .main-footer {
                margin-left: 250px; /* Align header and footer with the content */
            }
        }
      "))
    ),
    # Loading Screen
    div(
      id = "loading-screen",
      style = "height: 100vh; display: flex; justify-content: center; align-items: center; position: fixed; width: 100%; background-color: white; z-index: 1000;",
      div(
        class = "spinner-border",
        style = "width: 3rem; height: 3rem; border-width: 0.25em;",
        role = "status"
      ),
      h2("Loading...", style = "margin-left: 10px;")
    ),
    
    div(textOutput("rt_message")),
    div(textOutput("error_message")),
    
    # Actual app content
    div(
      id = "app-content",
      style = "visibility: hidden;",
      

        tabItems(
          tabItem(
            tabName = "files",
            
            fluidRow(
              bs4Card(
                title = "Feature table file",
                collapsible = TRUE,
                closable = FALSE,
                maximizable = TRUE,
                width = 12,
                
                div(
                  id = "FT_file",
                  #fileInput(inputId = "file1", "Choose a feature table file (.csv)"),
                  uiOutput("FTfile"),
                  uiOutput("picker"),
                  fluidRow(
                    column(6, checkboxInput("showTable", "Show feature table", value = FALSE)),
                    column(6, checkboxInput("showDemoTable", "Show demo feature table", value = FALSE))
                  ),
                  tags$hr(),
                  conditionalPanel(
                    condition = "input.showDemoTable == true",
                    tags$table(
                      class = "caption-top table",
                      style = "border-collapse: collapse; width: 100%;",
                      tags$thead(
                        tags$tr(
                          tags$th("File", style = "border: none;"),
                          tags$th("File Description", style = "border: none;"),
                          tags$th("Download", style = "border: none;")
                        )
                      ),
                      tags$tbody(
                        tags$tr(
                          tags$td(radioButtons("demoFTFile", label = NULL, choices = c("PA14 feature table" = "PA14example1"), selected = "PA14example1"), style = "border: none;"),
                          tags$td(
                            div(
                              class = "fixed-width",
                              tags$span("This is a LC-ESI-MS example dataset from an untargeted metabolomics experiment 
                                    investigating Pseudomonas aeruginosa perturbed with different antibiotics. The feature table was generated using XCMS. "
                                        #tags$a(href = "https://doi.org/10.1128/mSystems.00610-21", "Franke, R. et al.")
                              )
                            ),
                            style = "border: none;",
                          ),
                          
                          tags$td(downloadButton("download1", ""), style = "border: none;")
                        ),
                        tags$tr(
                          tags$td(radioButtons("demoFTFile", label = NULL, choices = c("Si11 feature table" = "Si11example2"), selected = NULL), style = "border: none;"),
                          tags$td(
                            div(
                              class = "fixed-width",
                              tags$span("This LC-ESI-MS dataset was generated from a mixture of 11 pure standards that are prone to producing in source fragments.
                                    The feature table was generated using MZmine 4."#,
                                        #tags$a(href = "https://pubs.acs.org/doi/10.1021/ac202450g", "Neumann, S. et al.")
                              )
                              
                            ),
                            style = "border: none;",
                          ),
                          
                          tags$td(downloadButton("download2", ""), style = "border: none;")
                        )
                      )
                    ),
                    fluidRow(
                      column(6, checkboxInput("viewDemoTable", "View demo feature table", value = FALSE)),
                      column(6, actionBttn("uploadDemoTable", "Select demo", value = FALSE))
                      
                    )
                    
                  )
                ),
                conditionalPanel(
                  condition = "input.showTable == true",
                  dataTableOutput("featureTable")
                ),
                conditionalPanel(
                  condition = "input.viewDemoTable == true",
                  dataTableOutput("demoFeatureTable")
                )
              )
            ),
            
            fluidRow(
              bs4Card(
                title ="MS2 file",
                collapsible = TRUE,
                closable = FALSE,
                maximizable = TRUE,
                width = 12,
                div(id = "MS2_file_menu",
                    shinycssloaders::withSpinner(uiOutput("MS2file")),# "Choose a pool MS2 file (.mzXML) or from MZmine .mgf") ,  
                    fluidRow(
                      column(6, checkboxInput("showMS2Table", "Show MS2 data table", value = FALSE)),
                      column(6, checkboxInput("showDemoMS2Table", "Show demo MS2 data table", value = FALSE))
                    ),
                    tags$hr(),
                    conditionalPanel(
                      condition = "input.showDemoMS2Table== true",
                      tags$table(
                        class = "caption-top table",
                        style = "border-collapse: collapse; width: 100%;",
                        tags$thead(
                          tags$tr(
                            tags$th("MS2 File", style = "border: none;"),
                            tags$th("MS2 File Description", style = "border: none;"),
                            tags$th("Download", style = "border: none;")
                          )
                        ),
                        tags$tbody(
                          tags$tr(
                            tags$td(radioButtons("demoMS2File", label = NULL, choices = c("PA14 MS2 file" = "MS2example1"), selected = NULL), style = "border: none;"),
                            tags$td(
                              div(
                                class = "fixed-width",
                                tags$span("DDA MS2 data in mzXML-format of a pool sample corresponding to the PA14 feature table."#,
                                          #tags$a(href = "https://doi.org/10.1128/mSystems.00610-21", "Franke, R. et al.")
                                )
                              ),
                              style = "border: none;",
                            ),
                            
                            tags$td(downloadButton("download3", ""), style = "border: none;")
                          ),
                          tags$tr(
                            tags$td(radioButtons("demoMS2File", label = NULL, choices = c("Si11 MS2 file" = "MS2example2"), selected = NULL), style = "border: none;"),
                            tags$td(
                              div(
                                class = "fixed-width",
                                tags$span("MGF file generated from DDA MS2 data of the Si11 mixture using MZmine 4."
                                )
                                
                              ),
                              style = "border: none;",
                            ),
                            
                            tags$td(downloadButton("download4", ""), style = "border: none;")
                          )
                        )
                      ),
                      fluidRow(
                        column(6, checkboxInput("viewMS2DemoTable", "View demo MS2 table", value = FALSE)),
                        column(6, actionBttn("UploadMS2DemoTable", "Select demo", value = FALSE)))
                    )
                ), #div
                
                
                conditionalPanel(
                  condition = "input.showMS2Table== true",
                  dataTableOutput("MS2table")),
                conditionalPanel(
                  condition = "input.viewMS2DemoTable== true",
                  dataTableOutput("demoMS2Table"))
              )
            ),
            
            fluidRow(
              bs4Card(title ="Library file",
                      collapsible = TRUE,
                      closable = FALSE,
                      maximizable = TRUE,
                      width = 12,
                      div(id = "Library_file_menu", 
                          
                          uiOutput("metabolites_file"),
                          
                          fluidRow(
                            column(6, checkboxInput("showMetabo", "Show metabolite target list", value = FALSE)),
                            column(6, checkboxInput("showDemoMetabo", "Show demo metabolite target list", value = FALSE))
                          ),
                          tags$hr(),
                          conditionalPanel(
                            condition = "input.showDemoMetabo== true",
                            tags$table(
                              class = "caption-top table",
                              style = "border-collapse: collapse; width: 100%;",
                              tags$thead(
                                tags$tr(
                                  tags$th("Metabolite Target List", style = "border: none;"),
                                  tags$th("File Description", style = "border: none;"),
                                  tags$th("Download", style = "border: none;")
                                )
                              ),
                              tags$tbody(
                                tags$tr(
                                  tags$td(radioButtons("demoMetaboFile", label = NULL, choices = c("PA14 Target List" = "PA14_target_example1"), selected = NULL), style = "border: none;"),
                                  tags$td(
                                    div(
                                      class = "fixed-width",
                                      tags$span("PA14 target list: A target list compiled from the Pseudomonas aeruginosa Metabolome Database (PAMDB)."
                                                #tags$a(href = "https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkx1061", "Wilks, A. et al.")
                                      )
                                    ),
                                    style = "border: none;",
                                  ),
                                  
                                  tags$td(downloadButton("download5", ""), style = "border: none;")
                                ),
                                tags$tr(
                                  tags$td(radioButtons("demoMetaboFile", label = NULL, choices = c("Si11 target file" = "Si11_target_example2"), selected = NULL), style = "border: none;"),
                                  tags$td(
                                    div(
                                      class = "fixed-width",
                                      tags$span("A target list of 11 standards."
                                                #tags$a(href = "https://pubs.acs.org/doi/10.1021/ac202450g", "Neumann, S. et al.")
                                      )
                                      
                                    ),
                                    style = "border: none;",
                                  ),
                                  
                                  tags$td(downloadButton("download6", ""), style = "border: none;")
                                )
                              )
                            ),
                            fluidRow(
                              column(6, checkboxInput("viewDemoMetabo", "Show metabolite target list", value = FALSE)),
                              column(6, actionBttn("UploadDemoMetabo", "Select demo", value = FALSE)))
                          )
                      ), #div
                      
                      conditionalPanel(
                        condition = "input.showMetabo== true",
                        dataTableOutput("metabo")),
                      conditionalPanel(
                        condition = "input.viewDemoMetabo== true",
                        dataTableOutput("demoMetaboTable"))
                      
              )
            ),
            
            fluidRow(
              bs4Card(title ="Neutral loss file",
                      collapsible = TRUE,
                      closable = FALSE,
                      maximizable = TRUE,
                      width = 12,
                      div(id = "NL_file_menu",
                          uiOutput("NL_file"),
                          checkboxInput("showNL", "Show Data Table", value = FALSE), 
                          DT::dataTableOutput("NLtable"))
                      
              )
            ),
            fluidRow(
              bs4Card(title ="ESI MS Adducts",
                      collapsible = TRUE,
                      closable = FALSE,
                      maximizable = TRUE,
                      width = 12,
                      div(id = "MS_Adducts",
                          uiOutput("adduct_file") ,
                          checkboxInput("showAdducts", "Show Data Table", value = FALSE), 
                          DT::dataTableOutput("adductsTable"))
                      
              )
            )
          ),
          
          # Tab 2: select input parameters
          tabItem(
            tabName = "parameters",
            fluidRow(
              box( title ="Filter the feature table",
                   collapsed=TRUE,
                   # slider bar: rt range ----
                   numericInput("minValue", "Minimum RT value (in seconds):", value = 60, min = 1, max = 1600),
                   numericInput("maxValue", "Maximum RT value (in seconds):", value = 1200, min = 1, max = 1600),
                   sliderInput("RTrange", "RT Range in seconds:",
                               min = 1, max = 1600,
                               value = c(60, 1200))
              ),
              
              box( title ="Check isotopes and multiple charge states", 
                   collapsed=TRUE,
                   checkboxInput(inputId="Isocheck", labe="Check C13 isotopes and multiple charge states", 
                                 value = TRUE)),
              box(title = "Correlation method", 
                  collapsed=TRUE,
                  selectInput(
                    inputId = "cor_method",
                    label = "Correlation method:",
                    input.cor_methods,
                    multiple = FALSE,
                    selectize = TRUE)
              ),
              #tags$hr(),
              box(title = "Correlation threshold",
                  collapsed=TRUE,
                  numericInput(
                    inputId = "cor_thr",
                    label = "Correlation threshold:",
                    value = 0.8,min = 0,
                    max = 1)
              ),
              box(title = "Retention time threshold (in second)",
                  collapsed=TRUE,
                  numericInput(
                    inputId = "rt_thr",
                    label = "RT threshold for correlation:",
                    value = 2,min = 0,
                    max = 30),
                  numericInput(
                    inputId = "rt_thr_exact",
                    label = "RT threshold for metabolites identification:",
                    value = 30,min = 0,
                    max = 1200),
                  numericInput(
                    inputId = "rt_thr_precursor",
                    label = "RT threshold for precursor feature matching:",
                    value = 20,min = 0,
                    max = 60),
                  numericInput(
                    inputId = "rt_thr_MS2",
                    label = "RT threshold for MS2 feature matching:",
                    value = 2,min = 0,
                    max = 60),
                  numericInput(
                    inputId = "rt_thr_NL",
                    label = "RT threshold for NL feature matching:",
                    value = 2,min = 0,
                    max = 60),
                  numericInput(
                    inputId = "rt_thr_adducts",
                    label = "RT threshold for adducts feature matching:",
                    value = 2,min = 0,
                    max = 10)
              ),
              box(title = "Ion polarity and primary ions",
                  collapsed=TRUE,
                  selectInput(
                    inputId = "IonPolarity",
                    label = "Ion polarity:",
                    c("pos","neg"),
                    selected=c("pos"),
                    multiple = FALSE,
                    selectize = TRUE),
                  selectInput(
                    inputId = "PIon",
                    label = "Primary ion:",
                    c("[M+H]+","[M+Na]+","[M-H]-"), 
                    selected=c("[M+H]+","[M+Na]+"),
                    multiple = TRUE,
                    selectize = TRUE)),
              box(title = "ppm",
                  collapsed=TRUE,
                  numericInput(
                    inputId = "ppm_exact",
                    label = "ppm for exact m/z matching:",
                    value = 5,min = 0,
                    max = 20),
                  
                  numericInput(
                    inputId = "ppm_precursor",
                    label = "ppm for precursor m/z matching:",
                    value = 5,min = 0,
                    max = 20),
                  
                  numericInput(
                    inputId = "ppm_MS2",
                    label = "ppm for MS2 m/z matching:",
                    value = 10,min = 0,
                    max = 20),
                  numericInput(
                    inputId = "ppm_NL",
                    label = "ppm for neutral losses m/z matching:",
                    value = 5,min = 0,
                    max = 20),
                  numericInput(
                    inputId = "ppm_adducts",
                    label = "ppm for adducts m/z matching:",
                    value = 5,min = 0,
                    max = 20)
                  
              ),
              box(title = "m/z difference tolerance", collapsed=TRUE,
                  numericInput(
                    inputId = "mz_diff_exact",
                    label = "m/z difference tolerance for exact m/z matching::",
                    value = 0.01,min = 0,
                    max = 1),
                  numericInput(
                    inputId = "mz_diff_precursor",
                    label = "m/z difference tolerance for precursor m/z matching:",
                    value = 0.01,min = 0,
                    max = 1),
                  numericInput(
                    inputId = "mz_diff_MS2",
                    label = "m/z difference tolerance for MS2 m/z matching:",
                    value = 0.01,min = 0,
                    max = 1),
                  
                  numericInput(
                    inputId = "mz_diff_NL",
                    label = "m/z difference tolerance for neutral losses m/z matching:",
                    value = 0.01,min = 0,
                    max = 1),
                  numericInput(
                    inputId = "mz_diff_adducts",
                    label = "m/z difference tolerance for adductsm/z matching:",
                    value = 0.01,min = 0,
                    max = 1)
                  
              ), 
              
              div(
                class = "col-12 text-right",  
                actionBttn("run", "Run",color = "primary")
              )
            ) #fluidRow
          ), # tab 2
          
          
          # Tab 3: output table
          tabItem(
            tabName = "output",
            selectInput("combinedFilter", "Filter by group or corgroup", choices = list("All" = "All")),
            fluidRow(
              bs4Card(
                title = "Output feature table",
                collapsible = TRUE,
                closable = FALSE,
                maximizable = TRUE,
                width = 12,
                id = "outputFeatureTableCard",
                shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "Output_FT")),
                div(
                  class = "row",
                  div(
                    class = "col-12 text-right",  # Adjust alignment to the right
                    downloadButton("downloadData", "")
                  )
                )
              )
            ),
            
            tags$head(
              tags$script(HTML(
                "document.addEventListener('DOMContentLoaded', function() {
      function setupResizeObserver(targetId) {
        var targetNode = document.getElementById(targetId);
        if (!targetNode) return;

        var observer = new MutationObserver(function(mutationsList) {
          for (var mutation of mutationsList) {
            if (mutation.type == 'attributes' && mutation.attributeName === 'class') {
              var isMaximized = $(mutation.target).hasClass('maximized-card');
              Shiny.setInputValue(targetId + 'Maximized', isMaximized, {priority: 'event'});

              if (isMaximized) {
                syncCardHeights();  // Adjust heights when maximized
              } else {
                resetCardHeights(); // Reset heights when minimized
              }
            }
          }
        });

        observer.observe(targetNode, { attributes: true });
      }

      function syncCardHeights() {
        var card1 = document.getElementById('networkPlotCard');
        var card2 = document.getElementById('boxPlotCard');

        if (card1 && card2) {
          var maxHeight = Math.max(card1.offsetHeight, card2.offsetHeight);
          card1.style.height = maxHeight + 'px';
          card2.style.height = maxHeight + 'px';
        }
      }

      function resetCardHeights() {
        var card1 = document.getElementById('networkPlotCard');
        var card2 = document.getElementById('boxPlotCard');

        if (card1 && card2) {
          card1.style.height = 'auto';
          card2.style.height = 'auto';
        }
      }

      // Setup resize observers for both bs4Cards
      setupResizeObserver('networkPlotCard');
      setupResizeObserver('boxPlotCard');
    });
    "
              ))
            ),

            tags$head(
              tags$style(HTML("
        .equal-height-row {
          display: flex;
          flex-wrap: wrap;
        }
        .equal-height-row .card {
          flex: 1;
          display: flex;
          flex-direction: column;
        }
        .equal-height-row .card-body {
          flex: 1;
        }
      "))
            ),
            fluidRow(
                     column(
                       width = 6,
                       bs4Card(
                         id = "networkPlotCard",
                         title = "Network plot",
                         collapsed = TRUE,
                         collapsible = TRUE,
                         closable = FALSE,
                         maximizable = TRUE,
                         width = NULL,  # Full width of the column
                         shinycssloaders::withSpinner(uiOutput("output_network")),
                         div(
                           style = "display:inline-block;margin-left: 92.5%",
                           downloadButton("downloadNetworkPlot", "")
                         )
                       )
                     ),
                     column(
                       width = 6,
                       bs4Card(
                         id = "boxPlotCard",
                         title = "Box plot",
                         collapsed = TRUE,
                         collapsible = TRUE,
                         closable = FALSE,
                         maximizable = TRUE,
                         width = NULL,  # Full width of the column
                         shinycssloaders::withSpinner(uiOutput("boxplot")),
                         div(
                           style = "display:inline-block;margin-left: 92.5%",
                           downloadButton("downloadboxPlot", "")
                         )
                       )
                     )
            ) ,# fluidRow
            fluidRow(
              bs4Card(
                id = "memoryUsageCard",
                title = "Memory usage",
                textOutput("memoryUsage"),
                width = 6
              )
            )# fluidRow
          ) # tab 3             
        )   #tabItems
    ) #div "app-content"
  ),#dashboardBody
  controlbar = NULL,  # Set controlbar to NULL
  footer = NULL,  # Set footer to NULL to avoid the footer error
  
  tags$head(
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        setTimeout(function(){
          $('#loading-screen').fadeOut('slow', function() {
            $('#app-content').css('visibility', 'visible');
          });
        }, 2000);  // Adjust the delay as needed
      });
    ")),
    tags$style(HTML("
      .spinner-border {
        display: inline-block;
        width: 2rem;
        height: 2rem;
        vertical-align: text-bottom;
        border: 0.25em solid currentColor;
        border-right-color: transparent;
        border-radius: 50%;
        animation: spinner-border 0.75s linear infinite;
      }

      @keyframes spinner-border {
        100% {
          transform: rotate(360deg);
        }
      }
    "))

  )
) # UI
