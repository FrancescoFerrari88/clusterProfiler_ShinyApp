# Shiny App for GO enrichment analysis with cluster profiler


Rlib = '/data/manke/group/ferrari/ShinyApps/Rlib_3.5'

.libPaths(Rlib)

library(shiny, lib.loc = Rlib)
library(ggplot2, lib.loc = Rlib)
library(clusterProfiler, lib.loc = Rlib)
library(VennDiagram, lib.loc = Rlib)
library(shinycssloaders, lib.loc = Rlib)
library(data.table, lib.loc = Rlib)
library(dplyr, lib.loc = Rlib)
require(org.Hs.eg.db, lib.loc = Rlib)
require(org.Mm.eg.db, lib.loc = Rlib)
require(org.Dm.eg.db, lib.loc = Rlib)

#library(shiny)
#library(ggplot2)
#library(clusterProfiler)
#library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
#library(shinycssloaders)
#library(data.table)
#library(dplyr)
#require(org.Hs.eg.db)
#require(org.Mm.eg.db)
#require(org.Dm.eg.db)

#################################
### USER INTERFACE DEFINITION ###
#################################

ui = navbarPage(strong("clusterProfiler GO Analysis "),
                
                # FIRST PANEL FOR DE VISUALIZATION                            
                tabPanel("Visualize DE Analysis",
                         
                         sidebarLayout(position="left",
                                       
                                       sidebarPanel(h3("Upload Panel"),
                                                    br(),
                                                    
                                                    fluidRow(
                                                      column(12,fileInput( "file", h4("Upload your Table"))),
                                                      br(),
                                                      br()
                                                      #column(3, checkboxInput("header", "header", TRUE))
                                                      #column(3, checkboxInput("index", "index", TRUE))
                                                    ),
                                                    
                                                    fluidRow(  
                                                      selectizeInput("select_ids", "Gene Names/IDs", 
                                                                     choices = NULL),
                                                      column(12, radioButtons("id_type", "ID type",
                                                                             choices = list("Symbol" = "SYMBOL", "Ensembl" = "ENSEMBL"
                                                                                            ),selected = "SYMBOL"))
                                                    ), 
                                                    
                                                    fluidRow(  
                                                      selectizeInput("select_expression", "Mean Expression", 
                                                                     choices = NULL)
                                                    ), 
                                                    
                                                    fluidRow(  
                                                      selectizeInput("select_logFC", "Log2 Fold Change", 
                                                                     choices = NULL)
                                                    ), 
                                                    
                                                    fluidRow(  
                                                      selectInput("select_pvalue", "pvalue", 
                                                                  choices = NULL)
                                                    ), 
                                                    
                                                    fluidRow(  
                                                      selectInput("select_fdr", "pvalue adjusted", 
                                                                  choices = NULL)
                                                    ), 
                                                    
                                                    fluidRow(  
                                                      sliderInput("fdr_thr", h5("pvalue adjusted threshold"),
                                                                  min = 0.0001, max = 0.1, value = 0.05)
                                                    ),
                                                    
                                                    fluidRow(  
                                                      sliderInput("fc_thr", h5("fold change threshold"),
                                                                  min = 1, max = 5, step=0.1, value = 1)
                                                    ),
                                                    
                                                    width = 3, align = "center"),
                                       
                                       mainPanel(
                                         tabsetPanel(
                                           tabPanel("Table", br(), dataTableOutput("Table")),
                                           
                                           tabPanel("MAplot - Volcano Plot",
                                                    br(),
                                                    sidebarLayout(position = "right",
                            
                                                                  sidebarPanel(h3("Summary"), br(), htmlOutput("summary"), width = 4, align = "center"),
                                                                  mainPanel(plotOutput("MAplot") %>% withSpinner(color="#0dc5c1"), plotOutput("Volcano") %>% withSpinner(color="#0dc5c1"))
                                                    )
                                           )
                                         )
                                       ))),
                
                
                
                
                # SECOND PANEL FOR FUNCTIONAL ENRICHMENT ANALYSIS
                tabPanel("Standard Enrichment Analysis",
                         
                         tabsetPanel(
                           tabPanel("Input Data ",
                                    
                                    br(), 
                                    
                                    sidebarLayout(position="right",
                                                  
                                                  sidebarPanel(h3("Control Panel"), align="center", width=3,
                                                               br(),
                                                               tabsetPanel(
                                                                 tabPanel("Input Data ",
                                                                          
                                                                          br(), 
                                                                          
                                                                          fluidRow(
                                                                            selectInput("model_organism", "Model Organism", 
                                                                                        choices = list("Homo sapiens" = "org.Hs.eg.db", 
                                                                                                       "Mus musculus" = "org.Mm.eg.db",
                                                                                                       "Drosophila melanogaster" = "org.Dm.eg.db"),
                                                                                        selected = "org.Mm.eg.db")
                                                                          ),
                                                                          
                                                                         # fluidRow(
                                                                         #   fileInput( "gtf", "Upload your GTF of reference (optional)")
                                                                         #  ),
                                                                          
    
                                                                          fluidRow(
                                                                            selectInput("target", "Target Set", 
                                                                                        choices = list("All DE Genes" = "All", 
                                                                                                       "Upregulated Genes" = "up",
                                                                                                       "Downregulated Genes" = "down"),
                                                                                        selected = "All")
                                                                          ),
                                                                          
                                                                          fluidRow(
                                                                            selectInput("background", "Background Set", 
                                                                                        choices = list("Expressed Genes" = "back_expr", 
                                                                                                       "Whole Genome" = "back_genome"),
                                                                                        selected = "back_expr")
                                                                          ),
                                                                          
                                                                          br()
                                                                          
                                                                 ),
                                                                 tabPanel("Analysis Parameters", 
                                                                          
                                                                          br(),
                                                                          
                                                                          fluidRow(
                                                                            selectInput("type_analysis", "Type of Analysis", 
                                                                                        choices = list("GO over-representation test" = "EnrichGO"
                                                                                                       #"GO classification" = "groupGO",
                                                                                                       #"KEGG over-represeantion test" = "EnrichKEGG",
                                                                                                       #"Disease Analysis" = "EnrichDO",
                                                                                                       #"Reactome Pathway Analysis" = "EnrichPathway",
                                                                                                       #"DAVID Functional Analysis" = "EnrichDAVID",
                                                                                                       #"Universal Enrichment Analisys" = "enricher"
                                                                                        ),
                                                                                        selected = "EnrichGO")
                                                                          ),
                                                                          
                                                                          fluidRow(
                                                                            selectInput("ont", "Subontology", 
                                                                                        choices = list("Molecular Functions" = "MF",
                                                                                                       "Biological Processes" = "BP", 
                                                                                                       "Cellular Compartments" = "CC",
                                                                                                       "All" = "All"
                                                                                        ),
                                                                                        selected = "BP")
                                                                          ),
                                                                          
                                                                          fluidRow(
                                                                            column(6, numericInput("pvalue_cutoff", "P-value threshold", value = 0.05)),
                                                                            column(6, numericInput("qvalue_cutoff", "Q-value threshold", value = 0.05))
                                                                            
                                                                          ),
                                                                          
                                                                          fluidRow(
                                                                            selectInput("pAdjustMethod", "P-value Adjust Method", 
                                                                                        choices = list("Benjamini-Hochberg" = "BH",
                                                                                                       "Bonferroni" = "bonferroni",
                                                                                                       "False Discovery Rate" = 'fdr',
                                                                                                       "Holm" = "holm",
                                                                                                       "Hochberg" = "hochberg",
                                                                                                       "Hommel" = "hommel",
                                                                                                       "Benjamini-Yekutieli"="BY"
                                                                                        ),
                                                                                        selected = "BH")
                                                                          )
                                                                          
                                                                 )
                                                                 
                                                               ),
                                                               
                                                               actionButton("submit_GO","Submit GO Analysis")
                                                               
                                                               
                                                  ),
                                                  
                                                  mainPanel(column(7, plotOutput("GO_input")), align="center"   #imageOutput("GO_input")
                                            
                                                            )
                                                  
                                                  )
                                    ),

                           tabPanel("Visualization",
                                    br(),
                                             navlistPanel(
          
                                               tabPanel("Barplot", 
                                                        sidebarLayout(position="right",
                                                                      sidebarPanel(align="center", width=3,
                                                                                  numericInput("num_show_1", h4("Number of GO terms to show"), value = 12)),
                                                        mainPanel(withSpinner(plotOutput("barplot"), color="#0dc5c1")))),
                                               tabPanel("Dotplot", 
                                                        sidebarLayout(position="right",
                                                                      sidebarPanel(align="center", width=3,
                                                                                   numericInput("num_show_2", h4("Number of GO terms to show"), value = 12)),
                                                        mainPanel(withSpinner(plotOutput("dotplot"), color="#0dc5c1")))),
                                               tabPanel("Emapplot", 
                                                        sidebarLayout(position="right",
                                                                      sidebarPanel(align="center", width=3,
                                                                                   numericInput("num_show_3", h4("Number of GO terms to show"), value = 12)),
                                                        mainPanel(withSpinner(plotOutput("emapplot"), color="#0dc5c1")))),
                                               
                                               #tabPanel("plotGOgraph", withSpinner(plotOutput("plotgograph"), color="#0dc5c1")),
                                              widths = c(2,8))      
                                    ),
                           tabPanel("Result Table", br(),downloadButton('downloadData', 'Download'),br(),br(), dataTableOutput("GO_output_table") %>% withSpinner(color="#0dc5c1")))),
 
               
                
                
                # THIRD PANEL FOR GSEA - TO BE CREATED
                #tabPanel("GSEA"
                #),
                
                
                # FOURTH PANEL FOR SESSIONINFO
                tabPanel("sessionInfo",
                         verbatimTextOutput("sessionInfo")
                )
)


#########################
### SERVER DEFINITION ###
#########################

server = function(input, output, session){
  data_in = reactive({
    req(input$file)
    #read.csv(input$file$datapath, sep="\t",header=input$header)
    as.data.frame(fread(input$file$datapath))
  })
  
  data_plot = reactive({
    req(input$file)
    #df = read.csv(input$file$datapath, sep="\t", header=input$header)
    df = as.data.frame(fread(input$file$datapath))
    df = df[,c(input$select_ids,input$select_expression, input$select_logFC, input$select_pvalue, input$select_fdr)]
    colnames(df) = c("GeneID","baseMean","log2FC","pvalue","padj")
        #handle gencode IDs
    if (input$id_type == "ENSEMBL"){
        if(grepl("\\.[0-9]{1,2}",df$GeneID[1])){df$GeneID<-gsub("\\.[0-9]+","",df$GeneID)}
      }
    df$Log_BaseMean = log10(df$baseMean)
    df$Log_pvalue = -log10(df$pvalue)
    df$sig = ifelse((df$padj < input$fdr_thr) & (abs(df$log2FC) > log2(input$fc_thr)), "Significant", "Not Significant")
    df$sig[is.na(df$padj)] = NA
    return(df)
  })
  
  output$Table = renderDataTable({
    data_in()
  })
  
  observeEvent(c(input$header, input$file), { 
    updateSelectizeInput(session, 'select_ids', choices = names(data_in()))
    updateSelectizeInput(session, 'select_expression', choices = names(data_in()))
    updateSelectizeInput(session, 'select_logFC', choices = names(data_in()))
    updateSelectizeInput(session, 'select_pvalue', choices = names(data_in()))
    updateSelectizeInput(session, 'select_fdr', choices = names(data_in()))
  })
  
  
  output$MAplot = renderPlot({
    ggplot(data_plot(), aes(Log_BaseMean, log2FC, color=sig)) + geom_point() + geom_hline(yintercept=0) + 
      geom_hline(yintercept=log2(input$fc_thr), linetype="dashed", color = "black") + 
      geom_hline(yintercept=-log2(input$fc_thr), linetype="dashed", color = "black")
  })
  
  output$Volcano = renderPlot({
    ggplot(data_plot(), aes(log2FC, Log_pvalue, color=sig)) + geom_point() + geom_vline(xintercept=0) +
      geom_vline(xintercept=log2(input$fc_thr), linetype="dashed", color = "black") + 
      geom_vline(xintercept=-log2(input$fc_thr), linetype="dashed", color = "black") 
  })
  
  output$summary = renderUI({
    df = na.omit(data_plot())
    str1 = paste0("Total significantly DE genes = ", sum(df$sig == "Significant"))
    str2 = paste0("Significantly  upregulated genes = ", sum(df$sig == "Significant" & df$log2FC > 0 ))       
    str3 = paste0("Significantly  downregulated genes = ", sum(df$sig == "Significant" & df$log2FC < 0 ))
    HTML(paste(str1, str2, str3, sep = '<br/><br/>'))
  })
  
  target_back = reactive({
    
    df = data_plot()
    
    ### define target gene set
    if (input$target == "All"){
      target = as.vector(df$GeneID[df$sig == "Significant"])
    }
    else if (input$target == "up"){
      target = as.vector(df$GeneID[df$sig == "Significant" & df$log2FC > 0])
    }
    else if (input$target == "down"){
      target = as.vector(df$GeneID[df$sig == "Significant" & df$log2FC < 0])
    }
    
    ### define background gene set
    if (input$background == "back_expr") {
      background_genes = as.vector(df$GeneID[df$sig == "Significant" | df$sig == "Not Significant"])
    }
    else {
      background_genes = as.vector(df$GeneID)
    }
    
    lista = list(na.omit(target),na.omit(background_genes))
    
    return(lista)
    
  })
  
  
  output$GO_input = renderPlot({      #renderImage({
    
    vd = venn.diagram(target_back(), 
                      category.names = c("target","backgr."), 
                      fill = c(1,4),
                      alpha = 0.3, 
                      filename = NULL,
                      lwd = 1,
                      lty = 'blank',
                      output = F ,
                      cex = 2,
                      cat.cex = 2,
                      resolution = 120)
    grid.newpage()
    grid.draw(vd)
  })
  
  
  
  ego_result = eventReactive(input$submit_GO, {
    
    #showNotification("Analysis started!")
    
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)
    require(org.Dm.eg.db)

    ### translate across IDs
    eg = bitr(target_back()[[1]], fromType=input$id_type, toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=input$model_organism)
    universe = bitr(target_back()[[2]], fromType=input$id_type, toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=input$model_organism)
    
    if (input$type_analysis == "EnrichGO"){
      
    ### run enrichGO 
      ego = enrichGO(
        gene = eg[["ENTREZID"]],
        universe = universe[["ENTREZID"]],
        OrgDb = input$model_organism[1],
        ont = input$ont,
        pAdjustMethod = input$pAdjustMethod,
        qvalueCutoff = input$qvalue_cutoff,
        pvalueCutoff = input$pvalue_cutoff,
        readable = TRUE )
      
    ### return enrichGO object for visualization
    return(ego)
    }
    
    #showNotification("Analysis is completed!")
    
  })  
  
  
  output$barplot = renderPlot({
      barplot(ego_result(), drop=T, showCategory=input$num_show_1)
  })
  
  output$dotplot = renderPlot({
      dotplot(ego_result(), showCategory=input$num_show_2)
  })
  
  output$emapplot = renderPlot({
      emapplot(ego_result(), showCategory=input$num_show_3)
  })
  
  #output$plotgograph = renderPlot({
  #  plotGOgraph(ego_result())
  #})
  
  output$GO_output_table = renderDataTable({
    as.data.frame(ego_result())
  })
  
  output$sessionInfo <- renderText({
    paste(capture.output(sessionInfo()),collapse = "\n")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Results_GO-Analysis_",Sys.Date(),".tsv", sep="")
    },
    content = function(file) {
      write.csv(as.data.frame(ego_result()), file, quote = F, sep="\t")
    }
  )
  
  #session$onSessionEnded(function() {
  #  junk <- dir(path="./www", pattern=".*.log") 
  #  file.remove(paste0("./www/",junk))
  #})
  
}



#################
### SHINY APP ###
#################

shinyApp(ui=ui,server=server)

