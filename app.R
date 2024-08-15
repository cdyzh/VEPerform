library(shiny)
library(dplyr)
library(yogiroc)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)

# Define UI for the application
ui <- fluidPage(
  titlePanel("VEPerform"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("upload_guide", "Upload Dataset Guide", class = "btn-primary"),
      selectInput("upload_type", "Select Upload Option:",
                  choices = list("Full Dataset (All Columns)" = "full", "Gene and Variant ID Only" = "gene_variant"),
                  selected = "full"),
      conditionalPanel(
        condition = "input.upload_type == 'full'",
        fileInput("file_full", "Upload Full CSV File", 
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),
      conditionalPanel(
        condition = "input.upload_type == 'gene_variant'",
        fileInput("file_gene_variant", "Upload Gene and Variant ID CSV File", 
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),
      selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
      checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
      checkboxGroupInput("scores", "Select Scores to Include:",
                         choices = list(
                           "VARITY",
                           "REVEL",
                           "AlphaMissense"
                         ),
                         selected = c("VARITY", "REVEL", "AlphaMissense")
      ),
      actionButton("plotButton", "Generate PRC Plot"),
      downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
      downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata"),
      helpText("Example genes: LDLR, TTN (no VARITY), ACADVL, DNM2, MYH7, SCN5A, etc."),
      helpText("If error occurs, it means there is not enough data to generate the PRC. Try using fewer predictors."), 
      actionButton("helpButton", "Help (Q&A)", class = "btn-warning")
    ),
    
    mainPanel(
      plotOutput("prcPlot", width = "600px", height = "600px"),
      textOutput("errorText")
    )
  )
)

server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  
  standard_colnames <- c(
    "base__hugo", 
    "gnomad__af", 
    "varity_r__varity_r", 
    "alphamissense__am_pathogenicity", 
    "revel__score", 
    "classification"
  )
  
  prcdata <- reactive({
    if (input$upload_type == "full" && !is.null(input$file_full)) {
      req(input$file_full)
      df <- read.csv(input$file_full$datapath, stringsAsFactors = FALSE)
      colnames(df) <- standard_colnames
      df
    } else if (input$upload_type == "gene_variant" && !is.null(input$file_gene_variant)) {
      req(input$file_gene_variant)
      df <- read.csv(input$file_gene_variant$datapath, stringsAsFactors = FALSE)
      
      # Load the full stored dataset for matching
      full_df <- read.csv("preprocessed_id.csv", stringsAsFactors = FALSE)
      colnames(df) <- c("base__hugo", "base__achange")
      
      # Merge based on gene and variant ID
      df <- merge(df, full_df, by = c("base__hugo", "base__achange"))
      print(df)
      df
    } else {
      read.table("preprocessed_id.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
    }
  })
  
  observe({
    df <- prcdata()
    if (!is.null(df)) {
      gene_names <- unique(df$base__hugo)
      updateSelectizeInput(session, "gene", choices = gene_names, selected = gene_names[1], server = TRUE)
    }
  })
  
  observe({
    req(input$gene)
    df <- prcdata()
    gene_data <- df %>% filter(base__hugo == input$gene)
    
    available_scores <- c()
    if (any(!is.na(gene_data$varity_r__varity_r))) {
      available_scores <- c(available_scores, "VARITY")
    }
    if (any(!is.na(gene_data$alphamissense__am_pathogenicity))) {
      available_scores <- c(available_scores, "AlphaMissense")
    }
    if (any(!is.na(gene_data$revel__score))) {
      available_scores <- c(available_scores, "REVEL")
    }
    
    updateCheckboxGroupInput(session, "scores", choices = available_scores, selected = available_scores)
  })
  
  observeEvent(input$plotButton, {
    df <- prcdata()
    gene_s <- input$gene
    exclude_common_variants <- input$common_variant_filter
    selected_scores <- input$scores
    
    names(df)[names(df) == "varity_r__varity_r"] <- "VARITY"
    names(df)[names(df) == "alphamissense__am_pathogenicity"] <- "AlphaMissense"
    names(df)[names(df) == "revel__score"] <- "REVEL"
    names(df)[names(df) == "gnomad__af"] <- "gnomAD_AF"
    
    if (exclude_common_variants) {
      df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
    }
    
    df <- df[order(df$classification), ]
    
    prcfiltered <- df %>%
      filter(base__hugo == gene_s)
    
    B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    
    tryCatch({
      yrobj <- yr2(truth = prcfiltered[["classification"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))
      
      plot_data(list(
        yrobj = yrobj,
        lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
        col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
        gene_s = gene_s,
        selected_scores = selected_scores,
        B_org = B_org,
        P_org = P_org,
        prcfiltered = prcfiltered  # Save the filtered data for metadata
      ))
      
      output$errorText <- renderText("")
      
    }, error = function(e) {
      plot_data(NULL)
      output$errorText <- renderText("Not enough data")
    })
    
    output$prcPlot <- renderPlot({
      plot_info <- plot_data()
      if (!is.null(plot_info)) {
        tryCatch({
          draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
          abline(h = 90, lty = "dashed")
          legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
        }, error = function(e) {
          showModal(modalDialog(
            title = 'Error',
            'Not enough data',
            easyClose = TRUE,
            footer = NULL
          ))
        })
      }
    }, width = 600, height = 600, res = 72)
  })
  
  output$downloadPlotPNG <- downloadHandler(
    filename = function() {
      paste("PRC_plot_", input$gene, ".png", sep = "")
    },
    content = function(file) {
      plot_info <- plot_data()
      if (!is.null(plot_info)) {
        png(file, width = 6, height = 6, units = "in", res = 72)
        draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
        abline(h = 90, lty = "dashed")
        legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
        dev.off()
      }
    }
  )
  
  output$downloadPlotPDF <- downloadHandler(
    filename = function() {
      paste("PRC_Report_", input$gene, ".pdf", sep = "")
    },
    content = function(file) {
      plot_info <- plot_data()
      if (!is.null(plot_info)) {
        # Generate the PDF report
        rmarkdown::render(input = "report_template.Rmd",
                          output_file = file,
                          params = list(
                            gene_s = plot_info$gene_s,
                            selected_scores = plot_info$selected_scores,
                            B_org = plot_info$B_org,
                            P_org = plot_info$P_org,
                            prcfiltered = plot_info$prcfiltered
                          ),
                          envir = new.env(parent = globalenv()))
      }
    }
  )
  
  observeEvent(input$upload_guide, {
    showModal(modalDialog(
      title = "Upload Dataset Guide",
      HTML("<h4>Option 1: Full Dataset (All Columns)</h4>
            <p>Upload a CSV file with the following columns:</p>
            <ul>
              <li><b>gene:</b> Gene name</li>
              <li><b>gnomad_af:</b> gnomAD allele frequency</li>
              <li><b>varity_r:</b> VARITY score</li>
              <li><b>alphamissense:</b> AlphaMissense score</li>
              <li><b>revel:</b> REVEL score</li>
              <li><b>classification:</b> TRUE for Pathogenic, FALSE for Benign</li>
            </ul>
            <h4>Option 2: Gene and Variant ID Only</h4>
            <p>Upload a CSV file with the following columns:</p>
            <ul>
              <li><b>gene:</b> Gene name</li>
              <li><b>variant_id:</b> Variant ID</li>
            </ul>
            <p>The data will be matched with the stored dataset to retrieve corresponding scores.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$helpButton, {
    showModal(modalDialog(
      title = "Q&A",
      HTML("<p><b>What is balanced precision?</b></p>
            <p>Balanced precision is useful in situations where the class distribution is imbalanced. In other words, it is the precision that would have been expected had the proportion of positive examples been balanced (equal to 50%).</p>
           <p>Definition at https://doi.org/10.1016/j.ajhg.2021.08.012</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
