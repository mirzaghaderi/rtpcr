#if (!requireNamespace("rtpcr", quietly = TRUE)) install.packages("rtpcr")
#if (!requireNamespace("rtpcr", quietly = TRUE)) install.packages("multcompView")
#if (!requireNamespace("markdown", quietly = TRUE)) install.packages("markdown")

library(multcompView)
library(rtpcr)
library(shiny)
library(dplyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(lme4)
library(lmerTest)
library(markdown)

# UI
ui <- fluidPage(
  titlePanel('rtpcr: qPCR Data Analysis & Plotting'),
  
  br(),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "func_tabs",
                  tabPanel("Home", helpText("Welcome! Please select an analysis tool from these tabs.")),
                  tabPanel("ANOVA_DCt",
                           radioButtons("src_dc", "Data Source", choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_dc == 'user'", fileInput("file_dc", "Upload CSV", accept = ".csv")),
                           numericInput("numFactors_dc", "Number of factors", 1, min = 1),
                           numericInput("numRefGenes_dc", "Number of reference genes", 1, min = 1),
                           textInput("block_dc", "Block column (optional)", ""),
                           checkboxInput("set40_dc", "Set missing Ct to 40", FALSE),
                           selectInput("pAdj_dc", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           actionButton("run_dc", "Run ANOVA_DCt")
                  ),
                  tabPanel("ANOVA_DDCt",
                           radioButtons("src_ddct", "Data Source", choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_ddct == 'user'", fileInput("file_ddct", "Upload CSV", accept = ".csv")),
                           numericInput("numFactors_ddct", "Number of factors", 1, min = 1),
                           numericInput("numRefGenes_ddct", "Number of reference genes", 1, min = 1),
                           # CHANGED: Replaced mainFactor column index with specs string
                           textInput("specs_ddct", "Specs (e.g., Concentration or Concentration | Type)", ""),
                           textInput("block_ddct", "Block column (optional)", ""),
                           textInput("calibrator_ddct", "Calibrator Level (NULL or single level; needs character!)", ""),
                           selectInput("seType_ddct", "SE type", choices = c("single.group",  "two.group", "paired.group")),
                           textInput("model_ddct", "model (optional)", ""),
                           selectInput("pAdj_ddct", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_ddct", "Set missing Ct to 40", FALSE),
                           actionButton("run_ddct", "Run ANOVA_DDCt")
                  ),
                  tabPanel("Efficiency",
                           radioButtons("src_eff", "Data Source",choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_eff == 'user'", fileInput("file_eff", "Upload CSV", accept = ".csv")),
                           numericInput("baseSize_eff", "Base font size", 12, min = 6),
                           textInput("extra_eff", "Further adjustments (e.g. + theme_bw() + labs(...))", ""),
                           hr(),
                           h4("Save Settings"),
                           fluidRow(
                             column(6, numericInput("eff_w", "Width (in)", 8, min = 1)),
                             column(6, numericInput("eff_h", "Height (in)", 6, min = 1))
                           ),
                           actionButton("run_eff", "Run Efficiency")
                  ),
                  tabPanel("TTEST_DDCt",
                           radioButtons("src_tt", "Data Source", choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_tt == 'user'", fileInput("file_tt", "Upload CSV", accept = ".csv")),
                           numericInput("numRefGenes_tt", "Number of reference genes", 1),
                           textInput("factorLevels_tt", "Factor levels (comma separated)", ""),
                           checkboxInput("paired_tt", "Paired test", FALSE),
                           checkboxInput("equalVar_tt", "Equal variance", TRUE),
                           selectInput("pAdj_tt", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_tt", "Set missing Ct to 40", FALSE),
                           actionButton("run_tt", "Run TTEST_DDCt")
                  ),
                  tabPanel("WILCOX_DDCt",
                           radioButtons("src_wx", "Data Source", choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_wx == 'user'", fileInput("file_wx", "Upload CSV", accept = ".csv")),
                           numericInput("numRefGenes_wx", "Number of reference genes", 1),
                           textInput("factorLevels_wx", "Factor levels (comma separated)", ""),
                           checkboxInput("paired_wx", "Paired test", FALSE),
                           selectInput("pAdj_wx", "p-value adjustment", choices = c("none","holm","bonferroni","fdr")),
                           checkboxInput("set40_wx", "Set missing Ct to 40", FALSE),
                           actionButton("run_wx", "Run WILCOX_DDCt")
                  ),
                  tabPanel("plotFactor",
                           radioButtons("src_pf", "Data Source", choices = c("Upload CSV" = "user", "Sample Data" = "sample"), selected = "user"),
                           conditionalPanel("input.src_pf == 'user'", fileInput("file_pf", "Upload CSV for Plotting", accept = ".csv")),
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
                           textInput("extra_pf", "Further adjustments (e.g. + theme_classic() + labs(...))", ""),
                           hr(),
                           h4("Save Settings"),
                           fluidRow(
                             column(6, numericInput("pf_w", "Width (in)", 8, min = 1)),
                             column(6, numericInput("pf_h", "Height (in)", 6, min = 1))
                           ),
                           actionButton("run_pf", "Run plotFactor")
                  )
      )
    ),
    mainPanel(
      tabsetPanel(id = "output_tabs",
                  tabPanel("Introduction",
                           div(style = "padding: 15px;",
                               includeMarkdown("README.md")
                           )
                  ),
                  tabPanel("ANOVA_DCt",
                           selectInput("gene_dc", "Select gene:", choices = character(0)),
                           tabsetPanel(id = "sub_dc",
                                       tabPanel("Input data", tableOutput("preview_dc")),
                                       tabPanel("Relative Expression", br(), downloadButton("download_dc","Download"), br(), tableOutput("table_dc")),
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
                  tabPanel("ANOVA_DDCt",
                           selectInput("gene_ddct", "Select gene:", choices = character(0)),
                           tabsetPanel(id = "sub_ddct",
                                       tabPanel("Input data", tableOutput("preview_ddct")),
                                       tabPanel("Relative Expression", br(), downloadButton("download_ddct", "Download"), br(), tableOutput("table_ddct")),
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
                  tabPanel("Efficiency",
                           tabsetPanel(id = "sub_eff",
                                       tabPanel("Input data", tableOutput("preview_eff")),
                                       tabPanel("Table", br(), downloadButton("download_eff", "Download CSV"), br(), tableOutput("table_eff")),
                                       tabPanel("Plot",
                                                plotOutput("plot_eff"),
                                                hr(),
                                                downloadButton("download_eff_png", "Download Plot as PNG"),
                                                downloadButton("download_eff_pdf", "Download Plot as PDF"))
                           )
                  ),
                  tabPanel("TTEST_DDCt",
                           tabsetPanel(id = "sub_tt",
                                       tabPanel("Input Data", tableOutput("preview_tt")),
                                       tabPanel("Results", br(), downloadButton("download_tt", "Download CSV"), br(), tableOutput("table_tt"))
                           )
                  ),
                  tabPanel("WILCOX_DDCt",
                           tabsetPanel(id = "sub_wx",
                                       tabPanel("Input Data", tableOutput("preview_wx") ),
                                       tabPanel("Results", br(), downloadButton("download_wx", "Download CSV"), br(), tableOutput("table_wx"))
                           )
                  ),
                  tabPanel("plotFactor",
                           tabsetPanel(id = "sub_pf",
                                       tabPanel("Input Data", tableOutput("preview_pf")),
                                       tabPanel("Plot",
                                                plotOutput("plot_pf_main", height = "600px"),
                                                hr(),
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
  
  # Logic to update ANOVA_DCt fields when Sample Data is selected
  observeEvent(input$src_dc, {
    if (input$src_dc == "sample") {
      updateNumericInput(session, "numFactors_dc", value = 3)
      updateNumericInput(session, "numRefGenes_dc", value = 1)
    }
  })
  
  # Logic to update ANOVA_DDCt fields when Sample Data is selected
  observeEvent(input$src_ddct, {
    if (input$src_ddct == "sample") {
      updateNumericInput(session, "numFactors_ddct", value = 2)
      updateNumericInput(session, "numRefGenes_ddct", value = 3)
      # CHANGED: Sample specs string
      updateTextInput(session, "specs_ddct", value = "Concentration")
      updateTextInput(session, "block_ddct", value = "block")
    }
  })
  
  # Logic to update plotFactor fields when Sample Data is selected
  observeEvent(input$src_pf, {
    if (input$src_pf == "sample") {
      updateSelectInput(session, "pf_group", selected = "gene")
      updateSelectInput(session, "pf_facet", selected = "gene")
      updateSelectInput(session, "pf_letters", selected = "sig")
    }
  })
  
  get_data <- function(src_input, file_input, sample_path) {
    if (src_input == "sample") {
      return(read.csv(sample_path))
    } else {
      req(file_input)
      return(read.csv(file_input$datapath))
    }
  }
  
  data_dc   <- reactive({ get_data(input$src_dc, input$file_dc, "exp/data_3factor.csv") })
  data_ddct <- reactive({ get_data(input$src_ddct, input$file_ddct, "exp/data_2factorBlock3ref.csv") })
  data_eff  <- reactive({ get_data(input$src_eff, input$file_eff, "eff/data_efficiency1.csv") })
  data_tt   <- reactive({ get_data(input$src_tt, input$file_tt, "exp/data_ttest18genes.csv") })
  data_wx   <- reactive({ get_data(input$src_wx, input$file_wx, "exp/data_ttest18genes.csv") })
  df_pf     <- reactive({ get_data(input$src_pf, input$file_pf, "plot/out_ANOVA_DDCt.csv") })
  
  output$preview_dc <- renderTable({ req(data_dc()); head(data_dc(), 50) })
  output$preview_ddct <- renderTable({ req(data_ddct()); head(data_ddct(), 50) })
  output$preview_eff <- renderTable({ req(data_eff()); head(data_eff(), 50) })
  output$preview_tt <- renderTable({ req(data_tt()); head(data_tt(), 50) })
  output$preview_wx <- renderTable({ req(data_wx()); head(data_wx(), 50) })
  output$preview_pf <- renderTable({ req(df_pf()); head(df_pf(), 50) })
  
  observeEvent(input$func_tabs, {
    tab_mapping <- c("ANOVA_DCt" = "ANOVA_DCt", "ANOVA_DDCt" = "ANOVA_DDCt",
                     "Efficiency" = "Efficiency", "TTEST_DDCt" = "TTEST_DDCt",
                     "WILCOX_DDCt" = "WILCOX_DDCt", "plotFactor" = "plotFactor")
    
    if (input$func_tabs %in% names(tab_mapping)) {
      updateTabsetPanel(session, "output_tabs", selected = tab_mapping[[input$func_tabs]])
    } else if (input$func_tabs == "Home") {
      updateTabsetPanel(session, "output_tabs", selected = "Introduction")
    }
  })
  
  observeEvent(input$run_dc,   { updateTabsetPanel(session, "sub_dc", selected = "Relative Expression") })
  observeEvent(input$run_ddct, { updateTabsetPanel(session, "sub_ddct", selected = "Relative Expression") })
  observeEvent(input$run_eff,  { updateTabsetPanel(session, "sub_eff", selected = "Table") })
  observeEvent(input$run_tt,   { updateTabsetPanel(session, "sub_tt", selected = "Results") })
  observeEvent(input$run_wx,   { updateTabsetPanel(session, "sub_wx", selected = "Results") })
  observeEvent(input$run_pf,   { updateTabsetPanel(session, "sub_pf", selected = "Plot") })
  
  res_eff <- eventReactive(input$run_eff, {
    res <- efficiency(data_eff(), base_size = input$baseSize_eff)
    if (nzchar(trimws(input$extra_eff))) {
      tryCatch({
        res$plot <- eval(parse(text = paste0("res$plot", input$extra_eff)))
      }, error = function(e) {
        showNotification(paste("Plot adjustment error:", e$message), type = "error")
      })
    }
    res
  })
  output$table_eff <- renderTable({ res_eff()$Efficiency })
  output$plot_eff <- renderPlot({ res_eff()$plot })
  
  res_tt <- eventReactive(input$run_tt, {
    TTEST_DDCt(data_tt(), numberOfrefGenes = input$numRefGenes_tt,
               Factor.level.order = if(input$factorLevels_tt == "") NULL else trimws(unlist(strsplit(input$factorLevels_tt, ","))),
               paired = input$paired_tt, var.equal = input$equalVar_tt, p.adj = input$pAdj_tt, set_missing_target_Ct_to_40 = input$set40_tt)
  })
  output$table_tt <- renderTable({ res_tt() })
  
  res_wx <- eventReactive(input$run_wx, {
    tryCatch({ WILCOX_DDCt(data_wx(), numberOfrefGenes = input$numRefGenes_wx,
                           Factor.level.order = if (input$factorLevels_wx == "") NULL else trimws(unlist(strsplit(input$factorLevels_wx, ","))),
                           paired = input$paired_wx, p.adj = input$pAdj_wx, set_missing_target_Ct_to_40 = input$set40_wx)
    }, error = function(e) { showNotification(e$message, type = "error"); NULL })
  })
  output$table_wx <- renderTable({ res_wx() })
  
  observeEvent(df_pf(), {
    cols <- colnames(df_pf())
    pick_col <- function(i) cols[min(i, length(cols))]
    updateSelectInput(session, "pf_x", choices = cols, selected = pick_col(2))
    updateSelectInput(session, "pf_y", choices = cols, selected = pick_col(4))
    updateSelectInput(session, "pf_low", choices = cols, selected = pick_col(9))
    updateSelectInput(session, "pf_up", choices = cols, selected = pick_col(10))
    updateSelectInput(session, "pf_group", choices = c("None" = "", cols),
                      selected = if(input$src_pf == "sample") "gene" else "")
    updateSelectInput(session, "pf_facet", choices = c("None" = "", cols),
                      selected = if(input$src_pf == "sample") "gene" else "")
    updateSelectInput(session, "pf_letters", choices = c("None" = "", cols),
                      selected = if(input$src_pf == "sample") "sig" else "")
  })
  
  pf_plot_obj <- eventReactive(input$run_pf, {
    p <- plotFactor(data = df_pf(), x_col = input$pf_x, y_col = input$pf_y, Lower.se_col = input$pf_low,
                    Upper.se_col = input$pf_up, group_col = if(input$pf_group == "") NULL else input$pf_group,
                    facet_col = if(input$pf_facet == "") NULL else input$pf_facet,
                    letters_col = if(input$pf_letters == "") NULL else input$pf_letters,
                    letters_d = input$pf_letters_d, col_width = input$pf_col_width, err_width = input$pf_err_width,
                    dodge_width = input$pf_dodge_width, fill_colors = if (input$pf_fill_colors == "") NULL else trimws(unlist(strsplit(input$pf_fill_colors, ","))),
                    color = input$pf_color, alpha = input$pf_alpha, base_size = input$pf_base_size,
                    legend_position = if (input$pf_legend_none) "none" else c(input$pf_legend_x, input$pf_legend_y))
    
    if (nzchar(trimws(input$extra_pf))) {
      tryCatch({
        p <- eval(parse(text = paste0("p", input$extra_pf)))
      }, error = function(e) {
        showNotification(paste("Plot adjustment error:", e$message), type = "error")
      })
    }
    p
  })
  output$plot_pf_main <- renderPlot({ req(pf_plot_obj()); pf_plot_obj() })
  
  res_dc <- eventReactive(input$run_dc, {
    tryCatch({ ANOVA_DCt(data_dc(), numOfFactors = input$numFactors_dc, numberOfrefGenes = input$numRefGenes_dc,
                         block = if(input$block_dc == "") NULL else input$block_dc, p.adj = input$pAdj_dc,
                         set_missing_target_Ct_to_40 = input$set40_dc)
    }, error = function(e) { showNotification(e$message, type="error"); NULL })
  })
  observe({
    genes <- if(!is.null(res_dc()$perGene)) names(res_dc()$perGene) else character(0)
    updateSelectInput(session, "gene_dc", choices = genes)
  })
  output$table_dc <- renderTable({ res_dc()$relativeExpression })
  output$anova_dc <- renderPrint({ req(res_dc(), input$gene_dc); res_dc()$perGene[[input$gene_dc]]$ANOVA_table })
  output$formula_dc <- renderPrint({ req(res_dc(), input$gene_dc); res_dc()$perGene[[input$gene_dc]]$lm_formula })
  output$final_data_dc <- renderTable({ req(res_dc(), input$gene_dc); res_dc()$perGene[[input$gene_dc]]$Final_data })
  
  res_ddct <- eventReactive(input$run_ddct, {
    tryCatch({ ANOVA_DDCt(data_ddct(), numOfFactors = input$numFactors_ddct, numberOfrefGenes = input$numRefGenes_ddct,
                          specs = input$specs_ddct, block = if(input$block_ddct == "") NULL else input$block_ddct,
                          model = if(input$model_ddct == "") NULL else input$model_ddct,
                          calibratorLevel = if(input$calibrator_ddct == "") NULL else input$calibrator_ddct,
                          p.adj = input$pAdj_ddct, set_missing_target_Ct_to_40 = input$set40_ddct, se.type = input$seType_ddct)
    }, error = function(e) { showNotification(e$message, type="error"); NULL })
  })
  observe({
    genes <- if(!is.null(res_ddct()$perGene)) names(res_ddct()$perGene) else character(0)
    updateSelectInput(session, "gene_ddct", choices = genes)
  })
  output$table_ddct <- renderTable({ res_ddct()$relativeExpression })
  output$anova_ddct <- renderPrint({ req(res_ddct(), input$gene_ddct); res_ddct()$perGene[[input$gene_ddct]]$ANOVA_table })
  output$formula_ddct <- renderPrint({ req(res_ddct(), input$gene_ddct); res_ddct()$perGene[[input$gene_ddct]]$lm_formula })
  output$final_data_ddct <- renderTable({ req(res_ddct(), input$gene_ddct); res_ddct()$perGene[[input$gene_ddct]]$Final_data })
  
  output$download_dc <- downloadHandler(filename = "ANOVA_DCt.csv", content = function(f) write.csv(res_dc()$relativeExpression, f, row.names = FALSE))
  output$download_ddct <- downloadHandler(filename = "ANOVA_DDCt.csv", content = function(f) write.csv(res_ddct()$relativeExpression, f, row.names=F))
  output$download_eff <- downloadHandler(filename = "Efficiency.csv", content = function(f) write.csv(res_eff()$Efficiency, f, row.names=F))
  output$download_tt <- downloadHandler(filename = "TTEST.csv", content = function(f) write.csv(res_tt(), f, row.names=F))
  output$download_wx <- downloadHandler(filename = "WILCOX_DDCt.csv", content = function(f) write.csv(res_wx(), f, row.names = FALSE))
  
  output$download_eff_png <- downloadHandler(
    filename = "Efficiency_Plot.png",
    content = function(file) ggsave(file, plot = res_eff()$plot, device = "png", width = input$eff_w, height = input$eff_h)
  )
  output$download_eff_pdf <- downloadHandler(
    filename = "Efficiency_Plot.pdf",
    content = function(file) ggsave(file, plot = res_eff()$plot, device = "pdf", width = input$eff_w, height = input$eff_h)
  )
  
  output$download_pf_plot <- downloadHandler(
    filename = "qPCR_Plot.png",
    content = function(file) ggsave(file, plot = pf_plot_obj(), device = "png", width = input$pf_w, height = input$pf_h)
  )
  output$download_pf_pdf <- downloadHandler(
    filename = "qPCR_Plot.pdf",
    content = function(file) ggsave(file, plot = pf_plot_obj(), device = "pdf", width = input$pf_w, height = input$pf_h)
  )
  
  output$download_lm_dc <- downloadHandler(filename = function() paste0("LM_", input$gene_dc, ".rds"), content = function(file) saveRDS(res_dc()$perGene[[input$gene_dc]]$lm, file))
  output$download_lm_ddct <- downloadHandler(filename = function() paste0("LM_", input$gene_ddct, ".rds"), content = function(file) saveRDS(res_ddct()$perGene[[input$gene_ddct]]$lm, file))
}

shinyApp(ui, server)