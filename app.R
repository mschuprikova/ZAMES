library(shiny)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(bslib)

# ============================================================
# UI
# ============================================================

ui <- navbarPage(
  
  title = "ZAMES 1.0",
  id = "zames_tabs",
  
  tabPanel(
    title = "Home",
    value = "home",
    
    fluidRow(
      
      column(
        width = 6,
        
        h2("Welcome to ZAMES"),
        
        includeMarkdown("home_page.md"),
        
        actionButton(
          inputId = "go_to_upload",
          label = "Go to ZAMES",
          class = "btn-success"
        )
      ),
      
      column(
        width = 6,
        h3("Principle of Action"),
         img(
           src = "zames-1.png",
           width = "100%"
         )
      )
    )
  ),
  
  tabPanel(
    title = "Upload to ZAMES",
    value = "upload",
    
    sidebarLayout(
      
      sidebarPanel(
        
        h4("Input data"),
        
        fileInput(
          inputId = "control_files",
          label = "Select control files (one or more)",
          multiple = TRUE,
          accept = ".ab1"
        ),
        
        fileInput(
          inputId = "sample_files",
          label = "Select sample files (one or more)",
          multiple = TRUE,
          accept = ".ab1"
        ),
        
        fileInput(
          inputId = "reference_file",
          label = "Select reference file",
          multiple = FALSE,
          accept = c(".fasta", ".fna", ".fa")
        ),
        
        textInput(
          inputId = "target_region",
          label = "Input Target region sequence",
          value = "For example: TTGTGTGGAATTGTGAGCGGATAA..."
        ),
        
        numericInput(
          inputId = "target_site",
          label = "Input Target site",
          value = 1,
          min = 1
        ),
        
        br(),
        
        actionButton(
          inputId = "check_input",
          label = "Check input",
          class = "btn-secondary"
        ),
        
        br(),
        br(),
        
        uiOutput("check_message"),
        
        br(),
        
        actionButton(
          inputId = "start_zames",
          label = "Start ZAMES",
          class = "btn-primary",
          width = "100%"
        )
      ),
      
      mainPanel(
        
        h3("How to get started with ZAMES?"),
        
        includeMarkdown("start_zames.md"),
        
        #img(
        #  src = "zames_1.png",
        #  width = "80%"
        #),
        
        br()
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  source("ZAMES_functions.R")
  
  imported <-  reactiveVal(NULL)
  results <- reactiveVal(NULL)
  selected_position <- reactiveVal(NULL)
  results_tab_created <- reactiveVal(FALSE)

  observeEvent(input$go_to_upload, {
    
    updateNavbarPage(
      session = session,
      inputId = "zames_tabs",
      selected = "upload"
    )
    
  })
  
  observeEvent(input$check_input, {
    
    errors <- character(0)
    
    # Control files
    if (is.null(input$control_files)) {
      errors <- c(errors, "Control files are not selected.")
    }
    
    # Sample files
    if (is.null(input$sample_files)) {
      errors <- c(errors, "Sample files are not selected.")
    }
    
    # Reference
    if (is.null(input$reference_file)) {
      errors <- c(errors, "Reference file is not selected.")
    }
    
    # Target region
    if (is.null(input$target_region) || input$target_region == "") {
      errors <- c(errors, "Target region is not specified.")
    }
    
    # Target site
    if (is.null(input$target_site) || input$target_site < 1) {
      errors <- c(errors, "Target site is not specified correctly.")
    }
    
    output$check_message <- renderUI({
      
      if (length(errors) == 0) {
        
        div(
          class = "alert alert-success",
          "Input data are ready."
        )
        
      } else {
        
        div(
          class = "alert alert-danger",
          lapply(
            errors,
            function(x) {
              div(x)
            }
          )
        )
      }
      
    })
    
  })

  observeEvent(input$start_zames, {
    
    req(input$control_files)
    req(input$sample_files)
    req(input$reference_file)
    req(input$target_region)
    req(input$target_site)
    
    disable("start_zames")
    
    withProgress(
      
      message = "Running ZAMES:",
      value = 0,
      
      {
  
        # 1.IMPORT
        incProgress(0.1, detail = "Importing data...")
        imported_data <- dataImportRun(input)
        
        # 2.READ DIRECTION
        incProgress(0.1, detail = "Checking read direction...")
        imported_data <- directionRun(imported_data)
        
        # 3.READ CONVERSION
        incProgress(0.1, detail = "Read conversing...")
        imported_data <- signalConversionRun(direction_info = imported_data, visualization = TRUE)
        
        # 4.SIGNAL CORRECTION
        incProgress(0.3, detail = "Correcting signals...")
        corrected_signals <- signalCorrectionRun(conversed_info = imported_data, visualization = TRUE)
        
        # 5.BACKGROUND THRESHOLD
        incProgress(0.1, detail = "Background Threshold...")
        background_threshold <- backgroundThresholdRun(imported_data = imported_data, corrected_signals = corrected_signals)
        
        # 6.MUTATION FREQUENCY
        incProgress(0.1, detail = "Mutation Frequency...")
        mutation_frequency <- mutationFrequencyRegion(imported_data = imported_data, background_threshold = background_threshold)
        
        # 7.REPORT PREPARATION
        incProgress(0.1, detail = "Mutation Frequency...")
        target_plots <- reportPreparationRun(
          imported_data = imported_data, # from step 3  
          corrected_signals = corrected_signals, # from step 4
          background_threshold = background_threshold, # from step 5
          mutation_frequency = mutation_frequency # from step 6
        )
        
        # 8.SAVE RESULTS 
        incProgress(0.1, detail = "Preparing results...")
        imported(imported_data)
        
        results(list(
          target_plots = target_plots,
          mutation_frequency = mutation_frequency
        ))
        
        # 9.FINISH ZAMES
        incProgress(1, detail = "Done!")
        
      }
    )

    if (!"results" %in% names(input$zames_tabs)) {
      
      insertTab(
        inputId = "zames_tabs",
        tabPanel(
          title = "Results",
          value = "results",
          
          h3("Chromatograms and area fractions"),
          uiOutput("target_plots_ui"),
          
          tags$hr(),
          
          h3("Mutation Frequency"),
          
          fluidRow(
            column(
              width = 6,
              h4("Target Region"),
              DT::DTOutput("mutation_summary")
            ),
            column(
              width = 6,
              h4(textOutput("selected_position_title", inline = TRUE)),
              DT::DTOutput("mutation_details")
            )
          )
        ),
        target = "upload",
        position = "after",
        select = TRUE,
        session = session
      )
      
    } else {
      
      updateNavbarPage(
        session = session,
        inputId = "zames_tabs",
        selected = "results"
      )
      
    }
    
    enable("start_zames")
    
  })
    
  output$target_plots_ui <- renderUI({
    
    req(results())
      
    plots <- results()$target_plots
      
    make_plot_column <- function(group_name, output_prefix) {
        
      group_plots <- plots[[group_name]]
        
      if (is.null(group_plots) || length(group_plots) == 0) {
        return(
          column(
            width = 6,
            h4(group_name),
            p("No plots available.")
          )
        )
      }
        
      column(
        width = 6,
          
        h4(group_name),
          
        lapply(seq_along(group_plots), function(i) {
            
          read_name <- names(group_plots)[i]
            
          if (is.null(read_name) || is.na(read_name) ||
              read_name == "") {
            read_name <- paste("Read", i)
          }
            
          tagList(
            h5(read_name),
              
            plotOutput(
              outputId = paste0(output_prefix, "_", i),
              height = "650px",
              width = "100%"
            ),
              
            tags$hr()
          )
        })
      )
    }
      
    fluidRow(
      make_plot_column("Control", "control_target"),
      make_plot_column("Sample", "sample_target")
    )
  })
    
  observe({
      
    req(results())
      
    plots <- results()$target_plots
      
    for (group_name in c("Control", "Sample")) {
        
      group_plots <- plots[[group_name]]
        
      if (is.null(group_plots)) next
        
      output_prefix <- if (group_name == "Control") {
        "control_target"
      } else {
        "sample_target"
      }
        
      for (i in seq_along(group_plots)) {
          
        local({
            
          ii <- i
          prefix <- output_prefix
          current_group <- group_name
            
          output[[paste0(prefix, "_", ii)]] <- renderPlot({
              
            req(results())
              
            results()$target_plots[[current_group]][[ii]]$plot
              
          })
        })
      }
    }
  })
    
  # Build summary table from the list of mutation-frequency tables
  buildTargetRegionSummary <- function(mf) {
      
    req(length(mf) > 0)
      
    # Keep the original order and duplicate names
    nucleotides <- names(mf)
      
    if (is.null(nucleotides)) {
      nucleotides <- vapply(
        mf,
        function(tbl) {
          if (
            is.data.frame(tbl) &&
            "Mutation Status" %in% names(tbl) &&
            any(tbl[["Mutation Status"]] == "reference nucleotide",
                na.rm = TRUE)
          ) {
            tbl[["Nucleotide"]][
              which(tbl[["Mutation Status"]] == "reference nucleotide")[1]
            ]
          } else {
            NA_character_
          }
        },
        character(1)
      )
    }
      
    # Extract detected mutations from each table
    statuses <- vapply(
      mf,
      function(tbl) {
          
        if (
          !is.data.frame(tbl) ||
          !"Mutation Status" %in% names(tbl)
        ) {
          return("\u2013")
        }
          
        mutation_status <- trimws(
          as.character(tbl[["Mutation Status"]])
        )
          
        # Exclude noise, reference nucleotide and empty values
        excluded_statuses <- c(
          "background noise",
          "reference nucleotide",
          "-",
          "",
          "NA"
        )
          
        detected <- mutation_status[
          !is.na(mutation_status) &
            !tolower(mutation_status) %in%
            tolower(excluded_statuses)
          ]
          
        detected <- unique(detected)
          
        if (length(detected) == 0) {
          return("\u2013")
        }
          
        if (length(detected) == 1) {
          paste0("mutation is detected: ", detected)
        } else {
          paste0(
            "mutations are detected: ",
            paste(detected, collapse = ", ")
          )
        }
      },
      character(1)
    )
      
    data.frame(
      Position = seq_along(mf),
      Nucleotide = nucleotides,
      Status = statuses,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
    
  observeEvent(results(), {
    
    req(results())
      
    mf <- results()$mutation_frequency
    req(length(mf) > 0)

    target_site <- imported()$reference_target_info$
      target_info$target_site
      
    if (
      length(target_site) != 1 ||
      is.na(target_site) ||
      !is.numeric(target_site) ||
      target_site < 1 ||
      target_site > length(mf)
    ) {
      target_site <- input$target_site
    }
      
    if (
      length(target_site) != 1 ||
      is.na(target_site) ||
      target_site < 1 ||
      target_site > length(mf)
    ) {
      target_site <- 1L
    }
      
    selected_position(as.integer(target_site))
      
  }, ignoreInit = FALSE)
    
  output$mutation_summary <- DT::renderDT({
      
    req(results())
      
    mf <- results()$mutation_frequency
    summary_df <- buildTargetRegionSummary(mf)
      
    target_site <- imported()$reference_target_info$
      target_info$target_site
      
    if (
      length(target_site) != 1 ||
      is.na(target_site) ||
      target_site < 1 ||
      target_site > nrow(summary_df)
    ) {
      target_site <- input$target_site
    }
      
    if (
      length(target_site) != 1 ||
      is.na(target_site) ||
      target_site < 1 ||
      target_site > nrow(summary_df)
    ) {
      target_site <- NA_integer_
    }
      
    # Positions with mutations outside the target site
    orange_positions <- summary_df$Position[
      summary_df$Position != target_site &
        !is.na(summary_df$Status) &
        summary_df$Status != "–"
    ]
      
    DT::datatable(
      summary_df[, c("Position", "Nucleotide", "Status"), drop = FALSE],
      rownames = FALSE,
      selection = "single",
      options = list(
        pageLength = 25,
        lengthChange = FALSE,
        searching = FALSE,
        ordering = FALSE,
        info = FALSE,
        dom = "t",
        scrollX = TRUE,
        columnDefs = list(
          list(width = "55px", targets = 0)
        ),
        na = "NA",
        rowCallback = DT::JS(
          sprintf(
            "function(row, data) {
          var position = Number(data[0]);
          var status = String(data[2] || '');

          if (position === %s) {
            $(row).css({
              'background-color': '#D9EAD3',
              'font-weight': 'bold'
            });
          } else if (
            status !== '\\u2013' &&
            status !== '' &&
            status !== 'NA'
          ) {
            $(row).css('background-color', '#FCE4D6');
          }
        }",
            if (is.na(target_site)) "null" else as.character(target_site)
          )
        )
      )
    )
      
  })
    
  observeEvent(input$mutation_summary_rows_selected, {
      
    selected_row <- input$mutation_summary_rows_selected
      
    if (length(selected_row) == 1) {
      selected_position(as.integer(selected_row))
    }
      
  })
    
  output$mutation_details <- DT::renderDT({
      
    req(results())
    req(selected_position())
      
    mf <- results()$mutation_frequency
    index <- selected_position()
      
    req(index >= 1, index <= length(mf))
      
    tbl <- mf[[index]]
    req(is.data.frame(tbl))
      
    display_tbl <- tbl
    
    mf_col <- "Mutation Frequency"
    
    if (mf_col %in% names(display_tbl)) {
      display_tbl[[mf_col]] <- as.character(display_tbl[[mf_col]])
      display_tbl[[mf_col]][is.na(display_tbl[[mf_col]])] <- "NA"
    }
    
    DT::datatable(
      display_tbl,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        ordering = FALSE,
        info = FALSE,
        dom = "t",
        scrollX = TRUE
        )
      )
    })
    
  output$selected_position_title <- renderText({
    req(results())
    req(selected_position())
      
    paste("Position", selected_position())
  })
    
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(
  ui = ui,
  server = server
)

