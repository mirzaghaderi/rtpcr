if (!requireNamespace("rtpcr", quietly = TRUE)) install.packages("rtpcr")
library(rtpcr)
library(shiny)
library(dplyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(lme4)
library(lmerTest)

# UI
ui <- fluidPage(
  titlePanel(HTML('rtpcr: qPCR Data Analysis & Plotting (<a href="https://github.com/mirzaghaderi/rtpcr" target="_blank">Help</a>)')),
  br(),
  div(
    style = "font-size: 16px; margin-top: -10px; margin-bottom: 15px;",
    HTML(
      "Citation:<br>
       Mirzaghaderi, G. (2025). rtpcr: a package for statistical analysis and graphical
       presentation of qPCR data in R. 
       <a href='https://doi.org/10.7717/peerj.20185' target='_blank'>
       PeerJ 13:e20185</a>."
    )
  ),
  br(),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "func_tabs",
                  tabPanel("ANOVA_DCt",
                           fileInput("file_dc", "Upload CSV", accept = ".csv"),
                           numericInput("numFactors_dc", "Number of factors", 1, min = 1),
                           numericInput("numRefGenes_dc", "Number of reference genes", 1, min = 1),
                           textInput("block_dc", "Block column (optional)", ""),
                           #checkboxInput("analyseAll_dc", "Analyse all targets", TRUE),
                           checkboxInput("set40_dc", "Set missing Ct to 40", FALSE),
                           selectInput("pAdj_dc", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           textInput("model_dc", "model (optional)", ""),
                           #checkboxInput("modelSE_dc", "Model-based SE", TRUE),
                           actionButton("run_dc", "Run ANOVA_DCt")
                  ),
                  tabPanel("ANOVA_DDCt",
                           fileInput("file_ddct", "Upload CSV", accept = ".csv"),
                           numericInput("numFactors_ddct", "Number of factors", 1, min = 1),
                           numericInput("numRefGenes_ddct", "Number of reference genes", 1, min = 1),
                           numericInput("mainFactor_ddct", "Main factor column", 1, min = 1),
                           textInput("block_ddct", "Block column (optional)", ""),
                           textInput("factorLevels_ddct", "Main factor level order (comma separated)", ""),
                           selectInput("seType_ddct", "SE type", choices = c("paired.group", "two.group", "single.group")),
                           textInput("model_ddct", "model (optional)", ""),
                           #checkboxInput("modelSE_ddct", "Model-based SE", TRUE),
                           selectInput("pAdj_ddct", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_ddct", "Set missing Ct to 40", FALSE),
                           actionButton("run_ddct", "Run ANOVA_DDCt")
                  ),
                  tabPanel("Efficiency",
                           fileInput("file_eff", "Upload CSV", accept = ".csv"),
                           numericInput("baseSize_eff", "Base font size", 12, min = 6),
                           actionButton("run_eff", "Run Efficiency")
                  ),
                  tabPanel("TTEST_DDCt",
                           fileInput("file_tt", "Upload CSV", accept = ".csv"),
                           numericInput("numRefGenes_tt", "Number of reference genes", 1),
                           textInput("factorLevels_tt", "Factor levels (comma separated)", ""),
                           checkboxInput("paired_tt", "Paired test", FALSE),
                           checkboxInput("equalVar_tt", "Equal variance", TRUE),
                           selectInput("pAdj_tt", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_tt", "Set missing Ct to 40", FALSE),
                           actionButton("run_tt", "Run TTEST_DDCt")
                  ),
                  tabPanel("WILCOX_DDCt",
                           fileInput("file_wx", "Upload CSV", accept = ".csv"),
                           numericInput("numRefGenes_wx", "Number of reference genes", 1),
                           textInput("factorLevels_wx", "Factor levels (comma separated)", ""),
                           checkboxInput("paired_wx", "Paired test", FALSE),
                           selectInput("pAdj_wx", "p-value adjustment",
                                       choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_wx", "Set missing Ct to 40", FALSE),
                           actionButton("run_wx", "Run WILCOX_DDCt")
                  ),
                  tabPanel("plotFactor",
                           fileInput("file_pf", "Upload CSV for Plotting", accept = ".csv"),
                           selectInput("pf_x", "X Axis Column", choices = NULL),
                           selectInput("pf_y", "Y Axis Column", choices = NULL),
                           selectInput("pf_low", "Lower Error Column", choices = NULL),
                           selectInput("pf_up", "Upper Error Column", choices = NULL),
                           selectInput("pf_group", "Group (Fill) Column", choices = c("None" = "")),
                           selectInput("pf_facet", "Facet Column", choices = c("None" = "")),
                           selectInput("pf_letters", "Letters Column (Optional)", choices = c("None" = "")),
                           numericInput("pf_letters_d", "Letters vertical offset", 0.2, step = 0.05),
                           numericInput("pf_col_width", "Column width", 0.8, step = 0.05),
                           numericInput("pf_err_width", "Error bar width", 0.15, step = 0.05),
                           numericInput("pf_dodge_width", "Dodge width", 0.8, step = 0.05),
                           textInput("pf_fill_colors", "Fill colors (comma separated, optional)", ""),
                           textInput("pf_color", "Outline color", "black"),
                           sliderInput("pf_alpha", "Transparency (alpha)", min = 0, max = 1, value = 1, step = 0.05),
                           numericInput("pf_base_size", "Base font size", 12, min = 6),
                           checkboxInput("pf_legend_none", "Hide legend", FALSE),
                           fluidRow(column(6, numericInput("pf_legend_x", "Legend X position", 1, min = 0, max = 1, step = 0.05) ),
                             column(6, numericInput("pf_legend_y", "Legend Y position", 1,min = 0, max = 1, step = 0.05) ) ),
                           actionButton("run_pf", "Run plotFactor"),
                           
                  )
      )
    ),
    mainPanel(
      tabsetPanel(id = "output_tabs",
                  tabPanel("ANOVA_DCt Output",
                           selectInput("gene_dc", "Select gene:", choices = character(0)),
                           tabsetPanel(
                             tabPanel("Input data",
                                      tableOutput("preview_dc")
                             ),
                             tabPanel("Relative Expression",
                                      br(),
                                      downloadButton("download_dc","Download"),
                                      br(),
                                      tableOutput("table_dc")
                             ),
                             tabPanel("Per-Gene Results",
                                      tabsetPanel(
                                        tabPanel("ANOVA Table", verbatimTextOutput("anova_dc")),
                                        tabPanel("LM Object", downloadButton("download_lm_dc", "Download LM object (.rds)")),
                                        tabPanel("LM Formula", verbatimTextOutput("formula_dc")),
                                        tabPanel("Final table", tableOutput("final_data_dc"))
                                      )
                             )
                           )
                  ),
                  tabPanel("ANOVA_DDCt Output",
                           selectInput("gene_ddct", "Select gene:", choices = character(0)),
                           tabsetPanel(
                             tabPanel("Input data",
                                      tableOutput("preview_ddct")
                             ),
                             # Relative Expression tab with Download button on top
                             tabPanel("Relative Expression",
                                      br(),
                                      downloadButton("download_ddct", "Download"),
                                      br(),
                                      tableOutput("table_ddct")
                             ),
                             # Per-gene results remain unchanged
                             tabPanel("Per-Gene Results",
                                      tabsetPanel(
                                        tabPanel("ANOVA Table", verbatimTextOutput("anova_ddct")),
                                        tabPanel("LM Object", downloadButton("download_lm_ddct", "Download LM object (.rds)")),
                                        tabPanel("LM Formula", verbatimTextOutput("formula_ddct")),
                                        tabPanel("Final table", tableOutput("final_data_ddct"))
                                      )
                             )
                           )
                  ),
                  tabPanel("Efficiency Output",
                           tabsetPanel(
                             tabPanel("Input data",
                                      tableOutput("preview_eff")
                             ),
                             tabPanel(
                               "Table",
                               br(),
                               downloadButton("download_eff", "Download"),
                               br(),
                               tableOutput("table_eff")
                             ),
                             tabPanel("Plot", plotOutput("plot_eff"))
                           )
                  ),
                  tabPanel(
                    "TTEST_DDCt Output",
                    tabsetPanel(
                      tabPanel(
                        "Input Data",
                        tableOutput("preview_tt")
                      ),
                      tabPanel(
                        "Results",
                        br(),
                        downloadButton("download_tt", "Download CSV"),
                        br(),
                        tableOutput("table_tt")
                      )
                    )
                  ),
                  tabPanel("WILCOX_DDCt Output",
                    tabsetPanel(
                      tabPanel("Input Data", tableOutput("preview_wx") ),
                      tabPanel("Results",
                    br(),
                    downloadButton("download_wx", "Download CSV"),
                    br(),
                    tableOutput("table_wx")
                      )
                    )
                  ),
                  tabPanel("plotFactor Output",
                    tabsetPanel(
                      tabPanel("Input Data", tableOutput("preview_pf")),
                      tabPanel("Plot",
                        plotOutput("plot_pf_main", height = "600px"),
                        downloadButton("download_pf_plot", "Download Plot as PNG"),
                        downloadButton("download_pf_pdf", "Download Plot as PDF")
                      )
                    )
                  )
      )
    )
  )
)

# SERVER
server <- function(input, output, session) {

  output$preview_dc <- renderTable({
    req(data_dc())
    head(data_dc(), 50)   # optional: show first 50 rows only
  })
  output$preview_ddct <- renderTable({
    req(data_ddct())
    head(data_ddct(), 50)   # optional: show first 50 rows only
  })
  output$preview_eff <- renderTable({
    req(data_eff())
    head(data_eff(), 50)   # optional: show first 50 rows only
  })
  output$preview_tt <- renderTable({
    req(data_tt())
    head(data_tt(), 50)   # optional: show first 50 rows only
  })
  output$preview_wx <- renderTable({
    req(data_wx())
    head(data_wx(), 50)   # optional: show first 50 rows only
  })
  output$preview_pf <- renderTable({
    req(df_pf())
    head(df_pf(), 50)  # optional: first 50 rows
  })
  
  
  
  parse_colors <- function(x) {
    if (x == "") return(NULL) 
    trimws(unlist(strsplit(x, ",")))}
  get_legend_position <- reactive({
    if (input$pf_legend_none) {
      return("none")
    }
    c(input$pf_legend_x, input$pf_legend_y)
  })
  
  

  observeEvent(input$func_tabs, {
    # Map left sidebar tab names to the corresponding right output tab names
    tab_mapping <- c(
      "README" = "README Output",
      "ANOVA_DCt" = "ANOVA_DCt Output",
      "ANOVA_DDCt" = "ANOVA_DDCt Output",
      "Efficiency" = "Efficiency Output",
      "TTEST_DDCt" = "TTEST_DDCt Output",
      "WILCOX_DDCt" = "WILCOX_DDCt Output",
      "plotFactor" = "plotFactor Output"
    )
    # Get the corresponding right tab
    selected_tab <- tab_mapping[[input$func_tabs]]
    if (!is.null(selected_tab)) {
      updateTabsetPanel(session, "output_tabs", selected = selected_tab)
    }
  })


  trim_split <- function(x) { trimws(unlist(strsplit(x, ","))) }

  # plotFactor Logic
  df_pf <- reactive({ req(input$file_pf); read.csv(input$file_pf$datapath) })
  observeEvent(df_pf(), {
    cols <- colnames(df_pf())
    pick_col <- function(i) cols[min(i, length(cols))]
    updateSelectInput(session, "pf_x", choices = cols, selected = pick_col(2))
    updateSelectInput(session, "pf_y", choices = cols, selected = pick_col(4))
    updateSelectInput(session, "pf_low", choices = cols, selected = pick_col(9))
    updateSelectInput(session, "pf_up", choices = cols, selected = pick_col(10))
    updateSelectInput(session, "pf_group", choices = c("None" = "", cols), selected = "")
    updateSelectInput(session, "pf_facet", choices = c("None" = "", cols), selected = "")
    updateSelectInput(session, "pf_letters", choices = c("None" = "", cols), selected = "")
  })


  pf_plot_obj <- eventReactive(input$run_pf, {
    plotFactor(
      data = df_pf(),
      x_col = input$pf_x,
      y_col = input$pf_y,
      Lower.se_col = input$pf_low,
      Upper.se_col = input$pf_up,
      group_col = if(input$pf_group == "") NULL else input$pf_group,
      facet_col = if(input$pf_facet == "") NULL else input$pf_facet,
      letters_col = if(input$pf_letters == "") NULL else input$pf_letters,
      etters_d = input$pf_letters_d,
      col_width = input$pf_col_width,
      err_width = input$pf_err_width,
      dodge_width = input$pf_dodge_width,
      fill_colors = parse_colors(input$pf_fill_colors),
      color = input$pf_color,
      alpha = input$pf_alpha,
      base_size = input$pf_base_size,
      legend_position = get_legend_position()
    )
  })

  output$plot_pf_main <- renderPlot({ req(pf_plot_obj()); pf_plot_obj() })

  output$download_pf_plot <- downloadHandler(
    filename = function() "qPCR_Plot.png",
    content = function(file) { ggsave(file, plot = pf_plot_obj(), device = "png", width = 8, height = 6) }
  )
  output$download_pf_pdf <- downloadHandler(
    filename = function() "qPCR_Plot.pdf", content = function(file) { ggsave(file, plot = pf_plot_obj(), device = "pdf", width = 8, height = 6) }
  )


  output$download_lm_ddct <- downloadHandler(
    filename = function() {
      paste0("LM_", input$gene_ddct, ".rds")
    },
    content = function(file) {
      req(res_ddct(), input$gene_ddct)
      gene_obj <- res_ddct()$perGene[[input$gene_ddct]]
      if (is.null(gene_obj) || is.null(gene_obj$lm)) {
        stop("No LM object available")
      }
      saveRDS(gene_obj$lm, file)
    } )


  # ANOVA_DCt
  data_dc <- reactive({ req(input$file_dc); read.csv(input$file_dc$datapath) })
  
  res_dc <- eventReactive(input$run_dc, {
    tryCatch({
      ANOVA_DCt(
        data_dc(),
        numOfFactors = input$numFactors_dc,
        numberOfrefGenes = input$numRefGenes_dc,
        block = if(input$block_dc == "") NULL else input$block_dc,
        model = if(input$model_dc == "") NULL else input$model_dc,
        p.adj = input$pAdj_dc,
        #analyseAllTarget = input$analyseAll_dc,
        set_missing_target_Ct_to_40 = input$set40_dc,
        #modelBased_se = input$modelSE_dc
      )
    }, error = function(e) {
      showNotification(e$message, type="error")
      NULL
    })
  })
  observe({
    genes <- if(!is.null(res_dc()$perGene)) names(res_dc()$perGene) else character(0)
    updateSelectInput(session, "gene_dc", choices = genes)
  })
  output$table_dc <- renderTable({ res_dc()$relativeExpression })
  output$anova_dc <- renderPrint({
    req(res_dc(), input$gene_dc)
    res_dc()$perGene[[input$gene_dc]]$ANOVA_table
  })
  
  output$formula_dc <- renderPrint({
    req(res_dc(), input$gene_dc)
    res_dc()$perGene[[input$gene_dc]]$lm_formula
  })
  output$final_data_dc <- renderTable({
    req(res_dc(), input$gene_dc)
    res_dc()$perGene[[input$gene_dc]]$Final_data
  })
  output$download_dc <- downloadHandler(
    filename = "ANOVA_DCt.csv",
    content = function(f) {
      write.csv(res_dc()$relativeExpression, f, row.names = FALSE)
    }
  )
  output$download_lm_dc <- downloadHandler(
    filename = function() paste0("LM_", input$gene_dc, ".rds"),
    content = function(file) {
      req(res_dc(), input$gene_dc)
      gene_obj <- res_dc()$perGene[[input$gene_dc]]
      saveRDS(gene_obj$lm, file)
    } )
  

  # --- ANOVA_DDCt ---
  data_ddct <- reactive({ req(input$file_ddct); read.csv(input$file_ddct$datapath) })
  res_ddct <- eventReactive(input$run_ddct, {
    tryCatch({
      ANOVA_DDCt(data_ddct(), numOfFactors = input$numFactors_ddct, numberOfrefGenes = input$numRefGenes_ddct,
                 mainFactor.column = input$mainFactor_ddct,
                 block = if(input$block_ddct == "") NULL else input$block_ddct,
                 model = if(input$model_ddct == "") NULL else input$model_ddct,
                 mainFactor.level.order = if(input$factorLevels_ddct == "") NULL else trim_split(input$factorLevels_ddct),
                 p.adj = input$pAdj_ddct,
                 set_missing_target_Ct_to_40 = input$set40_ddct,
                 se.type = input$seType_ddct,
                 #modelBased_se = input$modelSE_ddct
                 )
    }, error = function(e) { showNotification(e$message, type="error"); NULL })
  })
  observe({
    genes <- if(!is.null(res_ddct()$perGene)) names(res_ddct()$perGene) else character(0)
    updateSelectInput(session, "gene_ddct", choices = genes)
  })
  output$table_ddct <- renderTable({ res_ddct()$relativeExpression })
  output$anova_ddct <- renderPrint({ res_ddct()$perGene[[input$gene_ddct]]$ANOVA_table })
  output$formula_ddct <- renderPrint({ res_ddct()$perGene[[input$gene_ddct]]$lm_formula })
  output$final_data_ddct <- renderTable({ res_ddct()$perGene[[input$gene_ddct]]$Final_data })
  output$download_ddct <- downloadHandler(filename = "ANOVA_DDCt.csv", content = function(f) write.csv(res_ddct()$relativeExpression, f, row.names=F))

  # Efficiency
  data_eff <- reactive({ req(input$file_eff); read.csv(input$file_eff$datapath) })
  res_eff <- eventReactive(input$run_eff, { efficiency(data_eff(), base_size = input$baseSize_eff) })
  output$table_eff <- renderTable({ res_eff()$Efficiency })
  output$plot_eff <- renderPlot({ res_eff()$plot })
  output$download_eff <- downloadHandler(filename = "Efficiency.csv", content = function(f) write.csv(res_eff()$Efficiency, f, row.names=F))

  # TTEST_DDCt
  data_tt <- reactive({ req(input$file_tt); read.csv(input$file_tt$datapath) })
  res_tt <- eventReactive(input$run_tt, {
    TTEST_DDCt(data_tt(), numberOfrefGenes = input$numRefGenes_tt,
               Factor.level.order = if(input$factorLevels_tt == "") NULL else trim_split(input$factorLevels_tt),
               paired = input$paired_tt, var.equal = input$equalVar_tt, p.adj = input$pAdj_tt, set_missing_target_Ct_to_40 = input$set40_tt)
  })
  output$table_tt <- renderTable({ res_tt() })
  output$download_tt <- downloadHandler(filename = "TTEST.csv", content = function(f) write.csv(res_tt(), f, row.names=F))

  
  
  data_wx <- reactive({
    req(input$file_wx)
    read.csv(input$file_wx$datapath)
  })
  
  res_wx <- eventReactive(input$run_wx, {
    tryCatch({
      WILCOX_DDCt(
        data_wx(),
        numberOfrefGenes = input$numRefGenes_wx,
        Factor.level.order =
          if (input$factorLevels_wx == "") NULL
        else trim_split(input$factorLevels_wx),
        paired = input$paired_wx,
        p.adj = input$pAdj_wx,
        set_missing_target_Ct_to_40 = input$set40_wx
      )
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
  })
  
  output$table_wx <- renderTable({
    res_wx()
  })
  
  output$download_wx <- downloadHandler(
    filename = "WILCOX_DDCt.csv",
    content = function(file) {
      write.csv(res_wx(), file, row.names = FALSE)
    }
  )
  
  
  
}

shinyApp(ui, server)
