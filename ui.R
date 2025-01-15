# ui.R
library(shiny)
library(shinyjs)
library(dplyr)
library(yogiroc)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)
library(httr)
library(jsonlite)
library(DT)

# Define UI for the application
ui <- fluidPage(
  
  useShinyjs(),  # Initialize shinyjs for resetting files after error
  # Title
  fluidRow(
    column(
      width = 12,
      align = "left",
      HTML("<h2 style='text-align:center; color:blue;'><strong>VEPerform</strong></h2>") # CHANGE title format
    )),
  
  # Navigation Bar
  navbarPage(
    title = NULL,
    id = "navbar",
    
    # Main VEPerform Tab
    tabPanel("Basic",
             fluidPage(
               # Introductory Text
               tags$h3("Welcome to VEPerform"),
               tags$p(HTML("
                      <p>How well do variant effect predictors (VEPs) work for your favorite gene?
                      Of course, performance depends on what VEP score threshold is used.</p>
                      <p>VEPerform helps you visualize this by plotting:</p>
                      <ul>
                        <li><strong>Precision</strong> (% of pathogenicity predictions that are correct) versus</li>
                        <li><strong>Recall</strong> (% of known pathogenic variants that were predicted)</li>
                      </ul>
                      <p>Select a gene to begin.</p>
                    ")),
               
               sidebarLayout(
                 sidebarPanel(
                   selectizeInput(
                     "Main_gene", "Gene of Interest:", choices = NULL, selected = character(0), options = list(maxOptions = 1000, placeholder = "Type gene name here...")
                   ),
                   
                   radioButtons("data_source", "Reference Set of Pathogenic/Benign:",
                                choices = list("Default (ClinVar unfiltered)" = "existing"),
                                selected = "existing"),
                   
                   helpText(HTML("
                    <p><strong>NOTE:</strong> Alternatively, upload your own reference set via the <strong style='color:blue;'>Advanced</strong> VEPerform tab.</p>")),
                   
                   checkboxInput("Main_common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                   checkboxGroupInput("Main_scores", "Select Scores to Include:",
                                      choices = list("VARITY", "REVEL", "AlphaMissense"),
                                      selected = c("VARITY", "REVEL", "AlphaMissense")),
                   actionButton("Main_plotButton", "Make Precision vs Recall Plot"),
                   uiOutput("Main_download_buttons")  # To group download buttons
                 ),
                 mainPanel(
                   plotOutput("Main_PRCPlot", width = "600px", height = "600px"),
                   uiOutput("thresholdTableUI"), # jumptag
                   textOutput("Main_ErrorText")
                 )
               )
             )
    ),
    
    # Advanced Tab
    tabPanel("Advanced", 
             fluidPage(
               # Introductory Text
               # tags$h3("Customizable PRC"),
               tags$p(HTML("
                  <p>With Advanced options you can:</p>
                  <ul>
                    <li>Upload a filtered CSV from Basic Mode</li>
                    <li>Upload your own reference variants with pathogenicity annotations</li>
                    <li>Add scores from additional VEPs:
                      <ul>
                        <li>Retrieved by VEPerform from OpenCRAVAT</li>
                        <li>Uploaded by you</li>
                      </ul>
                    </li>
                  </ul>
                ")),
               
               sidebarLayout(
                 sidebarPanel(
                   radioButtons("input_type", "Would you like to:",
                                choices = list("Use my own reference set / filtered set from Basic Mode" = "own",
                                               "Fetch VEP scores from OpenCRAVAT / add to my own reference set" = "fetch")),
                   
                   # For own file
                   conditionalPanel(
                     condition = "input.input_type == 'own'",
                     fileInput("file_full", "Upload Full Reference Set CSV File",
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                     actionButton("upload_guide", "Upload Guide", class = "btn-info")
                   ),
                   
                   # If fetch from OpenCRAVAT
                   conditionalPanel(
                     condition = "input.input_type == 'fetch'",
                     fileInput("file_fetch", "Upload CSV File. Ensure it contains columns chrom, pos, ref_base, and alt_base",
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                     downloadButton("download_template_oc", "Download CSV Template"),
                     helpText("Click 'Fetch VEP Data' after uploading your file to fetch data from OpenCRAVAT"),
                     
                     checkboxGroupInput("Fetch_scores", "Select what to fetch from OpenCRAVAT:",
                                        choices = c(
                                          "ClinVar" = "clinvar", 
                                          "GnomAD" = "gnomad", 
                                          "VARITY_R" = "varity_r", 
                                          "REVEL" = "revel", 
                                          "AlphaMissense" = "alphamissense", 
                                          "VARITY_ER" = "varity_er", 
                                          "PROVEAN" = "provean",
                                          "SIFT" = "sift",
                                          "PolyPhen-2" = "polyphen2"
                                        ),
                                        selected = c("clinvar", "gnomad")),
                     # Run OpenCRAVAT and fetch the scores
                     actionButton("fetchButton", "Fetch VEP Data")
                   ),
                   
                   # selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                   checkboxGroupInput("scores", "Select Scores to Include:", choices = NULL, selected = NULL),
                   # checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                   actionButton("plotButton", "Make Precision vs Recall Plot"),
                   uiOutput("download_buttons")  # To group download buttons
                 ),
                 
                 mainPanel(
                   plotOutput("prcPlot", width = "600px", height = "600px"),
                   textOutput("errorText")
                 )
               )
             )
             
    ),
    # About Tab
    tabPanel("About",
             fluidRow(
               column(12, offset = 0,
                      h2("About VEPerform"),
                      p("This tool allows you to evaluate the performance of variant effect predictors for your favorite genes."),
                      p("VEPerform was developed by Cindy Zhang at the Roth Lab at the University of Pittsburgh."),
                      p("Citation:")
               )
             )
    )
  )
)
