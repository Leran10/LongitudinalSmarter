#' Launch interactive model selection dashboard
#'
#' Opens a Shiny dashboard that helps users diagnose their data and
#' select appropriate statistical models for longitudinal analysis.
#'
#' @param data_matrix Optional: Pre-load a data matrix (features x samples)
#' @param metadata Optional: Pre-load metadata
#' @param launch.browser Logical; if TRUE, opens in browser (default: TRUE)
#'
#' @return Opens Shiny app (returns NULL invisibly when closed)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch empty dashboard to upload data
#' launch_model_selector()
#'
#' # Launch with pre-loaded data
#' launch_model_selector(data_matrix = my_counts, metadata = my_metadata)
#' }
launch_model_selector <- function(data_matrix = NULL, metadata = NULL, launch.browser = TRUE) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for the dashboard. Install with: install.packages('shiny')")
  }

  # UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("LongitudinalSmarter: Model Selection Dashboard"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::h4("Upload Data"),

        shiny::fileInput("count_file", "Upload Count/PA Matrix (CSV/TSV)",
                        accept = c(".csv", ".tsv", ".txt")),

        shiny::fileInput("metadata_file", "Upload Metadata (optional)",
                        accept = c(".csv", ".tsv", ".txt")),

        shiny::hr(),

        shiny::radioButtons("data_type", "Data Type",
                           choices = c("Auto-detect" = "auto",
                                     "Count/Abundance" = "count",
                                     "Presence/Absence" = "PA"),
                           selected = "auto"),

        shiny::hr(),

        shiny::actionButton("diagnose_btn", "Diagnose Data",
                           class = "btn-primary", width = "100%"),

        shiny::hr(),

        shiny::h5("Quick Guide:"),
        shiny::p("1. Upload your data matrix", style = "font-size: 12px;"),
        shiny::p("2. Click 'Diagnose Data'", style = "font-size: 12px;"),
        shiny::p("3. Review recommendations", style = "font-size: 12px;"),
        shiny::p("4. View diagnostic plots", style = "font-size: 12px;")
      ),

      shiny::mainPanel(
        width = 9,

        shiny::tabsetPanel(
          id = "main_tabs",

          # Welcome tab
          shiny::tabPanel(
            "Welcome",
            shiny::h3("Welcome to the LongitudinalSmarter Model Selector!"),
            shiny::br(),
            shiny::p("This dashboard helps you choose the right statistical model for your longitudinal microbiome data."),
            shiny::br(),
            shiny::h4("What does this tool do?"),
            shiny::tags$ul(
              shiny::tags$li("Analyzes your data characteristics (sparsity, overdispersion, zero-inflation)"),
              shiny::tags$li("Recommends appropriate statistical models"),
              shiny::tags$li("Generates diagnostic plots"),
              shiny::tags$li("Provides R code to run your analysis")
            ),
            shiny::br(),
            shiny::h4("How to use:"),
            shiny::tags$ol(
              shiny::tags$li("Upload your data matrix (features as rows, samples as columns)"),
              shiny::tags$li("Optionally upload metadata"),
              shiny::tags$li("Click 'Diagnose Data'"),
              shiny::tags$li("Review recommendations in the 'Recommendations' tab"),
              shiny::tags$li("Explore diagnostic plots in the 'Diagnostics' tab")
            ),
            shiny::br(),
            shiny::h4("Expected data format:"),
            shiny::p(shiny::strong("Data Matrix:")),
            shiny::tags$ul(
              shiny::tags$li("Rows = features (genes, contigs, OTUs, etc.)"),
              shiny::tags$li("Columns = samples"),
              shiny::tags$li("First column = feature IDs (optional row names)"),
              shiny::tags$li("Values = counts (for abundance) or 0/1 (for PA)")
            ),
            shiny::br(),
            shiny::p(shiny::strong("Metadata (optional):")),
            shiny::tags$ul(
              shiny::tags$li("Rows = samples (matching data matrix columns)"),
              shiny::tags$li("Columns = variables (diagnosis, timepoint, patient ID, etc.)")
            )
          ),

          # Data preview tab
          shiny::tabPanel(
            "Data Preview",
            shiny::h3("Uploaded Data Preview"),
            shiny::br(),
            shiny::h4("Data Matrix"),
            shiny::verbatimTextOutput("data_summary"),
            shiny::tableOutput("data_preview"),
            shiny::br(),
            shiny::h4("Metadata"),
            shiny::tableOutput("metadata_preview")
          ),

          # Recommendations tab
          shiny::tabPanel(
            "Recommendations",
            shiny::h3("Model Recommendations"),
            shiny::br(),
            shiny::verbatimTextOutput("recommendations_text"),
            shiny::br(),
            shiny::h4("Example R Code"),
            shiny::verbatimTextOutput("example_code")
          ),

          # Diagnostics tab
          shiny::tabPanel(
            "Diagnostics",
            shiny::h3("Data Diagnostics"),
            shiny::br(),
            shiny::plotOutput("diag_plot1", height = "400px"),
            shiny::br(),
            shiny::plotOutput("diag_plot2", height = "400px"),
            shiny::br(),
            shiny::plotOutput("diag_plot3", height = "600px")
          ),

          # Details tab
          shiny::tabPanel(
            "Detailed Metrics",
            shiny::h3("Detailed Diagnostic Metrics"),
            shiny::br(),
            shiny::verbatimTextOutput("detailed_diagnostics")
          )
        )
      )
    )
  )

  # Server
  server <- function(input, output, session) {

    # Reactive values
    rv <- shiny::reactiveValues(
      data = NULL,
      metadata = NULL,
      diagnosis = NULL
    )

    # Load pre-loaded data if provided
    shiny::observe({
      if (!is.null(data_matrix)) {
        rv$data <- data_matrix
      }
      if (!is.null(metadata)) {
        rv$metadata <- metadata
      }
    })

    # Read uploaded data
    shiny::observeEvent(input$count_file, {
      req(input$count_file)

      tryCatch({
        ext <- tools::file_ext(input$count_file$name)

        if (ext == "csv") {
          df <- read.csv(input$count_file$datapath, row.names = 1, check.names = FALSE)
        } else {
          df <- read.delim(input$count_file$datapath, row.names = 1, check.names = FALSE)
        }

        rv$data <- as.matrix(df)
        shiny::showNotification("Data uploaded successfully!", type = "message")

      }, error = function(e) {
        shiny::showNotification(paste("Error reading file:", e$message), type = "error")
      })
    })

    # Read metadata
    shiny::observeEvent(input$metadata_file, {
      req(input$metadata_file)

      tryCatch({
        ext <- tools::file_ext(input$metadata_file$name)

        if (ext == "csv") {
          df <- read.csv(input$metadata_file$datapath, row.names = 1, check.names = FALSE)
        } else {
          df <- read.delim(input$metadata_file$datapath, row.names = 1, check.names = FALSE)
        }

        rv$metadata <- df
        shiny::showNotification("Metadata uploaded successfully!", type = "message")

      }, error = function(e) {
        shiny::showNotification(paste("Error reading metadata:", e$message), type = "error")
      })
    })

    # Run diagnosis
    shiny::observeEvent(input$diagnose_btn, {
      req(rv$data)

      shiny::withProgress(message = "Analyzing data...", value = 0, {

        tryCatch({
          shiny::incProgress(0.3, detail = "Calculating statistics...")

          rv$diagnosis <- diagnose_data(
            data_matrix = rv$data,
            data_type = input$data_type,
            metadata = rv$metadata,
            make_plots = TRUE
          )

          shiny::incProgress(0.7, detail = "Generating plots...")

          shiny::showNotification("Diagnosis complete! Check the Recommendations tab.", type = "message")

          # Switch to recommendations tab
          shiny::updateTabsetPanel(session, "main_tabs", selected = "Recommendations")

        }, error = function(e) {
          shiny::showNotification(paste("Error during diagnosis:", e$message), type = "error")
        })
      })
    })

    # Data summary
    output$data_summary <- shiny::renderText({
      req(rv$data)

      paste0(
        "Dimensions: ", nrow(rv$data), " features × ", ncol(rv$data), " samples\n",
        "Value range: ", min(rv$data, na.rm = TRUE), " to ", max(rv$data, na.rm = TRUE), "\n",
        "Total observations: ", length(rv$data), "\n",
        "Missing values: ", sum(is.na(rv$data))
      )
    })

    # Data preview
    output$data_preview <- shiny::renderTable({
      req(rv$data)
      head(rv$data[, 1:min(6, ncol(rv$data)), drop = FALSE])
    }, rownames = TRUE)

    # Metadata preview
    output$metadata_preview <- shiny::renderTable({
      req(rv$metadata)
      head(rv$metadata)
    }, rownames = TRUE)

    # Recommendations
    output$recommendations_text <- shiny::renderPrint({
      req(rv$diagnosis)
      print(rv$diagnosis)
    })

    # Example code
    output$example_code <- shiny::renderText({
      req(rv$diagnosis)

      if (rv$diagnosis$data_type == "PA") {
        code <- paste0(
          "# Presence/Absence Analysis\n",
          "library(LongitudinalSmarter)\n\n",
          "results <- fit_longitudinal_PA(\n",
          "  PA_matrix = your_PA_matrix,\n",
          "  metadata = your_metadata,\n",
          "  dx_var = \"Dx.Status\",\n",
          "  timepoint_var = \"timepoint\",\n",
          "  random_var = \"patientID\",\n",
          "  covariates = c(\"Sex\", \"Age\")  # adjust as needed\n",
          ")\n\n",
          "# View results\n",
          "results$combined %>% filter(p.adj < 0.05)\n"
        )
      } else {
        fam <- rv$diagnosis$recommendations$recommended_families[1]
        if (grepl("PA", fam)) {
          code <- paste0(
            "# Data is very sparse - consider PA analysis\n",
            "library(LongitudinalSmarter)\n\n",
            "# Option 1: Convert to PA\n",
            "PA_matrix <- create_PA_matrix(\n",
            "  count_matrix = your_count_matrix,\n",
            "  target_cpm = 0.5,\n",
            "  min_reads_floor = 3\n",
            ")\n\n",
            "results <- fit_longitudinal_PA(\n",
            "  PA_matrix = PA_matrix,\n",
            "  metadata = your_metadata,\n",
            "  dx_var = \"Dx.Status\",\n",
            "  timepoint_var = \"timepoint\",\n",
            "  random_var = \"patientID\"\n",
            ")\n"
          )
        } else {
          code <- paste0(
            "# Count/Abundance Analysis\n",
            "library(LongitudinalSmarter)\n\n",
            "results <- fit_longitudinal_abundance(\n",
            "  count_matrix = your_count_matrix,\n",
            "  metadata = your_metadata,\n",
            "  dx_var = \"Dx.Status\",\n",
            "  timepoint_var = \"timepoint\",\n",
            "  random_var = \"patientID\",\n",
            "  covariates = c(\"Sex\", \"Age\"),  # adjust as needed\n",
            "  family = \"", fam, "\"  # Recommended based on your data\n",
            ")\n\n",
            "# View results\n",
            "results$combined %>% filter(p.adj < 0.05)\n"
          )
        }
      }

      code
    })

    # Diagnostic plots
    output$diag_plot1 <- shiny::renderPlot({
      req(rv$diagnosis)
      if (rv$diagnosis$data_type == "PA") {
        rv$diagnosis$plots$prevalence_distribution
      } else {
        rv$diagnosis$plots$zero_distribution
      }
    })

    output$diag_plot2 <- shiny::renderPlot({
      req(rv$diagnosis)
      if (rv$diagnosis$data_type == "count") {
        rv$diagnosis$plots$mean_variance
      }
    })

    output$diag_plot3 <- shiny::renderPlot({
      req(rv$diagnosis)
      if (rv$diagnosis$data_type == "count") {
        rv$diagnosis$plots$example_distributions
      }
    })

    # Detailed diagnostics
    output$detailed_diagnostics <- shiny::renderPrint({
      req(rv$diagnosis)
      str(rv$diagnosis$diagnostics, max.level = 2)
    })
  }

  # Run app
  shiny::shinyApp(ui = ui, server = server, options = list(launch.browser = launch.browser))
}
