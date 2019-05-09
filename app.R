
Rlib = '/data/manke/group/ferrari/ShinyApps/Rlib_3.5'
#Rlib = '/home/ferrari/R/x86_64-redhat-linux-gnu-library/3.5/'
# 
.libPaths(Rlib)
# 
library(shiny, lib.loc = Rlib)
library(shinyalert, lib.loc = Rlib)
library(shinydashboard, lib.loc = Rlib)
#library(dashboardthemes)
library(ggplot2, lib.loc = Rlib)
# library(clusterProfiler, lib.loc = Rlib)
library(VennDiagram, lib.loc = Rlib)
library(shinycssloaders, lib.loc = Rlib)
library(data.table, lib.loc = Rlib)
library(dplyr, lib.loc = Rlib)
#require(org.Hs.eg.db, lib.loc = Rlib)
#require(org.Mm.eg.db, lib.loc = Rlib)
#require(org.Dm.eg.db, lib.loc = Rlib)


futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")



ui <- dashboardPage(
  dashboardHeader(title="GO Analysis clusterProfiler", titleWidth = 300),
  dashboardSidebar(
    sidebarMenu(
      menuItem("  Load DE Analysis Table", tabName = "a", icon = icon("chart-line")),
      menuItem("  GO Enrichment Analysis", tabName = "b", icon = icon("dna")),
      menuItem("  Compare Clusters", tabName = "e", icon = icon("th-list")),
      menuItem("sessionInfo", tabName = "c", icon = icon("fingerprint")),
      menuItem("Documentation", tabName = "d", icon = icon("compass"))
    )
  ),
  dashboardBody(
    
    
    ### changing theme
    #shinyDashboardThemes(
    #  theme = "blue_gradient"
    #),
    
    tags$head(tags$style(HTML('
                              
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: white;
                              }
                                
                              

                              '))),
    
    tabItems(
      
      # First tab content
      tabItem(tabName = "a",
              
              sidebarLayout(position="left",
                            
                            box( title = "Upload Panel", status = "primary", solidHeader = TRUE,
                                 collapsible = TRUE, 
                                         br(),
                                 fileInput( "file", h4("Upload your Table")),
                                         br(),
                                 
                                 selectizeInput("select_ids", "Gene Names/IDs", 
                                                choices = NULL),
                                 radioButtons("id_type", "ID type",
                                                         choices = list("Symbol" = "SYMBOL", "Ensembl" = "ENSEMBL"
                                                         ),selected = "SYMBOL"),
                                 selectizeInput("select_expression", "Mean Expression", 
                                                choices = NULL),
                                         
                                 selectizeInput("select_logFC", "Log2 Fold Change", 
                                                choices = NULL),        
                            
                                 selectInput("select_pvalue", "pvalue", 
                                             choices = NULL),       
                                         
                                 selectInput("select_fdr", "pvalue adjusted", 
                                             choices = NULL),  
                                 
                                 useShinyalert(),
                                 actionButton("confirm_inp","Confirm Input Data"),
                                    
                                         
                                         width = 3, align = "center"),
                            
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Input Table", br(), dataTableOutput("Table")),
                                
                                tabPanel("MAplot - Volcano Plot",
                                         br(),
                                         sidebarLayout(position = "right",
                                                       
                                                       box(title = "Thresholds and Summary", status = "primary", solidHeader = TRUE,
                                                             collapsible = TRUE, 
            
                                                                    sliderInput("fdr_thr", h5("pvalue adjusted threshold"),
                                                                                min = 0.0001, max = 0.1, value = 0.05),
                                                                    
                                                                    sliderInput("fc_thr", h5("fold change threshold"),
                                                                                min = 1, max = 5, step=0.1, value = 1),
                                                           
                                                           br(), br(), htmlOutput("summary"), 
                                                           
                                                                    
                                                                    width = 4, align = "center"),
                                                       
                                                       mainPanel(plotOutput("MAplot") %>% withSpinner(color="#0dc5c1"), plotOutput("Volcano") %>% withSpinner(color="#0dc5c1"))
                                         )
                                ),
                                
                                
                                tabPanel("Selected Significant Genes", br(), dataTableOutput("Sel_Table")%>% withSpinner(color="#0dc5c1"))
                                
                              )
                            ))
      ),
      
      # Second tab content
      tabItem(tabName = "b",
              
              sidebarLayout(position="left",
                            
                            box(title = "Control Panel", status = "primary", solidHeader = TRUE,
                                                             collapsible = TRUE, align="center", width=3,
                                                                      br(),
                                                                      tabsetPanel(
                                                                        tabPanel("Input Data ",

                                                                                 br(),
                                                                                 selectInput("model_organism", "Model Organism",
                                                                                             choices = list("Homo sapiens" = "org.Hs.eg.db",
                                                                                                            "Mus musculus" = "org.Mm.eg.db",
                                                                                                            "Drosophila melanogaster" = "org.Dm.eg.db"),
                                                                                             selected = "org.Mm.eg.db"),


                                                                                 # fluidRow(
                                                                                 #   fileInput( "gtf", "Upload your GTF of reference (optional)")
                                                                                 #  ),

                                                                                 selectInput("target", "Target Set",
                                                                                             choices = list("All DE Genes" = "All",
                                                                                                            "Upregulated Genes" = "up",
                                                                                                            "Downregulated Genes" = "down"),
                                                                                             selected = "All"),

                                                                                 selectInput("background", "Background Set",
                                                                                             choices = list("Expressed Genes" = "back_expr",
                                                                                                            "All Genes" = "back_genome"),
                                                                                             selected = "back_expr"),

                                                                                 br()

                                                                        ),
                                                                        tabPanel("Parameters",

                                                                                 br(),
                                                                                 selectInput("type_analysis", "Type of Analysis",
                                                                                             choices = list("GO over-representation test" = "EnrichGO"
                                                                                                            #"GO classification" = "groupGO",
                                                                                                            #"KEGG over-represeantion test" = "EnrichKEGG",
                                                                                                            #"Disease Analysis" = "EnrichDO",
                                                                                                            #"Reactome Pathway Analysis" = "EnrichPathway",
                                                                                                            #"DAVID Functional Analysis" = "EnrichDAVID",
                                                                                                            #"Universal Enrichment Analisys" = "enricher"
                                                                                             ),
                                                                                             selected = "EnrichGO"),

                                                                                 selectInput("ont", "Subontology",
                                                                                             choices = list("Molecular Functions" = "MF",
                                                                                                            "Biological Processes" = "BP",
                                                                                                            "Cellular Compartments" = "CC",
                                                                                                            "All" = "All"
                                                                                             ),
                                                                                             selected = "BP"),

                                                                                 column(6, numericInput("pvalue_cutoff", "P-value threshold", value = 0.05)),
                                                                                 column(6, numericInput("qvalue_cutoff", "Q-value threshold", value = 0.05)),

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

                                                                      ),
                                                                      useShinyalert(),
                                                                      actionButton("submit_GO","Submit GO Analysis"),
                                                                      br(),
                                                                      br(),
                                                                      plotOutput("GO_input", height = "250px")
                                                                      #column(7, plotOutput("GO_input"))


                                                         ),
                            
                          
                            
                            mainPanel(
                              
                              tabsetPanel(
                                
                                tabPanel("Visualization",
                                                    br(),
                                                    navlistPanel(

                                                      tabPanel("Barplot",
                                                               sidebarLayout(position="right",

                                                                             box(title = "Number of GO terms to show", status = "primary", solidHeader = TRUE,
                                                                                 collapsible = F, align="center", width=4,
                                                                                 numericInput("num_show_1", "", value = 12)),

                                                                             mainPanel(width = 12, withSpinner(plotOutput("barplot"), color="#0dc5c1")))),
                                                      tabPanel("Dotplot",
                                                               sidebarLayout(position="right",

                                                                             box(title = "Number of GO terms to show", status = "primary", solidHeader = TRUE,
                                                                                 collapsible = F, align="center", width=4,
                                                                                 numericInput("num_show_2", "", value = 12)),

                                                                             mainPanel(width = 12, withSpinner(plotOutput("dotplot"), color="#0dc5c1")))),
                                                      tabPanel("Emapplot",
                                                               sidebarLayout(position="right",

                                                                             box(title = "Number of GO terms to show", status = "primary", solidHeader = TRUE,
                                                                                 collapsible = F, align="center", width=4,
                                                                                 numericInput("num_show_3", "", value = 12)),

                                                                             mainPanel(width = 12, withSpinner(plotOutput("emapplot"), color="#0dc5c1")))),


                                                      #tabPanel("plotGOgraph", withSpinner(plotOutput("plotgograph"), color="#0dc5c1")),
                                                      widths = c(2,8))
                                           ),
                                           tabPanel("Result Table", br(),
                                                    downloadButton('downloadData', 'Download'),
                                                    br(),
                                                    br(),
                                                    dataTableOutput("GO_output_table") %>% withSpinner(color="#0dc5c1"))

                              )
                              
                              
                            )
                            
                            )
              
      ),
      
      
      # Third tab content
      tabItem(tabName = "e",
              sidebarLayout(position="left",
                            
                            box(title = "Input Panel", status = "primary", solidHeader = TRUE,
                                collapsible = F, align="center", width=3,
                                br(),
                                selectInput("inp_Clust", "Import Clusters",
                                            choices = list("From DE Analysis (Up vs Down)" = 1, 
                                                           "From Table" = 2), 
                                                           #"Manual" = 3),
                                            selected = 1),
                                
                                
                                conditionalPanel(
                                  condition = "input.inp_Clust == 1",
                      
                                  selectInput("ont_cl", "Subontology", 
                                              choices = list("Molecular Functions" = "MF",
                                                             "Biological Processes" = "BP", 
                                                             "Cellular Compartments" = "CC",
                                                             "All" = "All"
                                              ),
                                              selected = "BP"),
                                  br(),
                                  useShinyalert(),
                                  actionButton("submit_GO_cl1","Submit Compare Cluster Analysis")
                                  ),
                                
                                
                                conditionalPanel(
                                  condition = "input.inp_Clust == 2",
                                  tags$hr(),
                                  fileInput( "file_clust", h4("Upload your Clustered Table")),
                                  selectizeInput("select_ids_cl", "Gene Names/IDs", 
                                                 choices = NULL),
                                  radioButtons("id_type_cl", "ID type",
                                               choices = list("Symbol" = "SYMBOL", "Ensembl" = "ENSEMBL"
                                               ),selected = "SYMBOL"),
                                  selectizeInput("clusters", "Clusters", 
                                                 choices = NULL),
                                  actionButton("confirm","Confirm Cluster Column"),
                                  br(),
                                  br(),

                                  checkboxGroupInput("clust_include", "Which clusters do you want to compare?",
                                                       choices=NULL),
                                  
                                  selectInput("model_organism_cl", "Model Organism", 
                                              choices = list("Homo sapiens" = "org.Hs.eg.db", 
                                                             "Mus musculus" = "org.Mm.eg.db",
                                                             "Drosophila melanogaster" = "org.Dm.eg.db"),
                                              selected = "org.Mm.eg.db"),
                                  selectInput("ont_cl", "Subontology", 
                                              choices = list("Molecular Functions" = "MF",
                                                             "Biological Processes" = "BP", 
                                                             "Cellular Compartments" = "CC",
                                                             "All" = "All"
                                              ),
                                              selected = "BP"),
                                  br(),
                                  actionButton("submit_GO_cl2","Submit Compare Cluster Analysis")
          
                                 ),
                                
                                conditionalPanel(
                                  condition = "input.inp_Clust == 3",
                                  tags$hr(),
                                  actionButton('addCluster', 'Add Cluster'),
                                  #actionButton('reset_cl', 'Reset'),
                                  br(),
                                  tags$div(id = 'placeholder'),
                                  selectInput("model_organism", "Model Organism", 
                                              choices = list("Homo sapiens" = "org.Hs.eg.db", 
                                                             "Mus musculus" = "org.Mm.eg.db",
                                                             "Drosophila melanogaster" = "org.Dm.eg.db"),
                                              selected = "org.Mm.eg.db"),
                                  selectInput("ont_cl", "Subontology", 
                                              choices = list("Molecular Functions" = "MF",
                                                             "Biological Processes" = "BP", 
                                                             "Cellular Compartments" = "CC",
                                                             "All" = "All"
                                              ),
                                              selected = "BP"),
                                  br(),
                                  actionButton("submit_GO_cl3","Submit Compare Cluster Analysis")
                                  
                                  
                                )
                                
                                
                                
                                ),
                            
                            mainPanel(width = 9,
                              
                              
                              
                              conditionalPanel(
                                condition = "input.inp_Clust == 1",
                                tabsetPanel(
                                  tabPanel("Visualize Results",
                                           br(),
                                           sidebarLayout(position="right",
                                                         
                                                         box(title = "Number of GO terms to show", status = "primary", solidHeader = TRUE,
                                                             collapsible = F, align="center", width=3,
                                                             numericInput("num_show_cl_DE", "", value = 12)),
                                                         
                                                         mainPanel(withSpinner(plotOutput("dotplot_cl_DE"), color="#0dc5c1")))
                                           
                                           ),
                                  
                                  tabPanel("Result Table", br(),
                                           downloadButton('downloadData_cl_DE', 'Download'),
                                           br(),br(), 
                                           dataTableOutput("output_table_cl_DE") %>% withSpinner(color="#0dc5c1"))
                                  
                                )),
                              
                              
                              conditionalPanel(
                                condition = "input.inp_Clust == 2",
                              tabsetPanel(
                                tabPanel("Input Table", 
                                         br(),
                                         dataTableOutput("cl_input")
                                         ),
                                tabPanel("Visualize Results", 
                                         br(),
                                           sidebarLayout(position="right",
                                                         
                                                         box(title = "Number of GO terms to show", status = "primary", solidHeader = TRUE,
                                                             collapsible = F, align="center", width=3,
                                                             numericInput("num_show_cl_Tab", "", value = 2)),
                                                         
                                                         mainPanel(withSpinner(plotOutput("comp_cl_Tab"), color="#0dc5c1"),width = "1000px")
                                )),
                                tabPanel("Result Table",
                                         br(),
                                         downloadButton('downloadData_cl_Tab', 'Download'),
                                         br(),
                                         br(),
                                         dataTableOutput("table_cl_Tab_result"))
                                
                                 )),

                            conditionalPanel(
                              condition = "input.inp_Clust == 3",
                              tabsetPanel(
                                tabPanel("Visualize Results"),
                                tabPanel("Result Table")
                              ))
                            
                            
                            )
      )),
     
      
      # Fourth tab content
      tabItem(tabName = "c",
              verbatimTextOutput("sessionInfo")
      ),
      
      
      # Fifth tab content
      tabItem(tabName = "d",
              tabPanel("Documentation",
                       br(),
                       downloadButton("DDS_v", 'Download Vignette'),
                       br(),
                       br(),
                       downloadButton("DDS", 'Download Sample Dataset (mouse)')
 
                       #tabPanel("Docs",
                       #        includeHTML("/data/manke/group/shiny/ferrari/clusterProfiler_GOenrich/ShinyApp_documentation/Documentation_ShinyApp_clusterProfiler.html")
                       #        )
                       
                       #htmlOutput("documentation") %>% withSpinner(color="#0dc5c1")
                       #includeHTML("/data/manke/group/shiny/ferrari/clusterProfiler_GOenrich/ShinyApp_documentation/Documentation_ShinyApp_clusterProfiler.html")
              )
      )
      
      
       
    )
  )
)

options(shiny.maxRequestSize=30*1024^2) 

server <- function(input, output, session) { 

  ###############################
  #### Visualize DE Analysis ####
  ###############################

  
  observeEvent(input$confirm_inp, {
    # Show a modal when the button is pressed
    shinyalert("PERFECT!", "Your data have been uploaded. Click on the 'MAplot-Volcano P  lot' tab to visualize your DE analysis and set significance thresholds. You can visualize and interrogate the set of significant genes by clicking on the 'Selected Significan Genes' tab!", type = "success")
  })
  
  
  data_in = reactive({
    req(input$file)
    library(clusterProfiler, lib.loc = Rlib)
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

  sel_tab = reactive({
    df_sig = data_plot()[data_plot()$sig=="Significant",c("GeneID","baseMean","log2FC","pvalue","padj")]
    df_sig = df_sig[rowSums(is.na(df_sig)) != ncol(df_sig),]
    return(df_sig)
  })
  output$Sel_Table = renderDataTable({
    sel_tab()
  })


  #observeEvent(c(input$header, input$file), {
  observeEvent(input$file, {
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


  ################################
  #### GO Enrichment Analysis ####
  ################################

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
                      height = 300, 
                      width = 300,
                      fill = c(1,4),
                      alpha = 0.3,
                      filename = NULL,
                      lwd = 1,
                      lty = 'blank',
                      output = F ,
                      cex = 1,
                      cat.cex = 1,
                      resolution = 120)
    grid.newpage()
    grid.draw(vd)
  })

  observeEvent(input$submit_GO, {
    # Show a modal when the button is pressed
    shinyalert("SUPER!", "The analysis will take few minutes to complete. The result will be visualized in this page. Click on the 'Result Table' tab to visualize and download the result table!", type = "success")
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

  #output$comp_up_down = renderPlot({
  #  if (!is.null(ego_result()[[2]])){
  #    dotplot(ego_result()[[2]], showCategory=input$num_show_3)
  #  }
  #})

  output$GO_output_table = renderDataTable({
    as.data.frame(ego_result())
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Results_GO-Analysis_",Sys.Date(),".tsv", sep="")
    },
    content = function(file) {
      write.csv(as.data.frame(ego_result()), file, quote = F, sep="\t")
    }
  )





  ##########################
  #### Compare Clusters ####
  ##########################

  observeEvent(input$submit_GO_cl1, {
    # Show a modal when the button is pressed
    shinyalert("GREAT!", "The analysis will take few minutes to complete. The result will be visualized in this page. Click on the 'Result Table' tab to visualize and download the result table!", type = "success")
  })
  
  observeEvent(input$submit_GO_cl2, {
    # Show a modal when the button is pressed
    shinyalert("AMAZING!", "The analysis will take few minutes to complete. Click on the 'Visualize Results' tab to visualize your results and on the 'Result Table' tab to visualize and download the result table!", type = "success")
  })
  
  data_in_clust = reactive({
    req(input$file_clust)
    library(clusterProfiler, lib.loc = Rlib)
    as.data.frame(fread(input$file_clust$datapath))
  })


  data_in_clust_analysis = reactive({
    req(input$file_clust)
    df_cl = as.data.frame(fread(input$file_clust$datapath))
    df_cl = df_cl[,c(input$select_ids_cl,input$clusters)]
    colnames(df_cl) = c("GeneID","Clusters")
    ### add handling of gencode ###
    if (input$id_type_cl == "ENSEMBL"){
      if(grepl("\\.[0-9]{1,2}",df_cl$GeneID[1])){df_cl$GeneID<-gsub("\\.[0-9]+","",df_cl$GeneID)}
    }

    return(df_cl)
  })


  output$cl_input = renderDataTable({
    data_in_clust_analysis()
  })


  clusts_EnrichGO = reactive({
    req(input$clusters)
    df_cl = data_in_clust_analysis()
    list_clust = list()
    for (i in 1:length(cl_final_analysis())){
      v = as.vector(df_cl$GeneID[df_cl$Clusters == cl_final_analysis()[i]])
      list_clust[[i]] = v
    }
    names(list_clust) = cl_final_analysis()
    return(list_clust)
  })


  observeEvent(input$file_clust, {
    updateSelectizeInput(session, 'select_ids_cl', choices = names(data_in_clust()))
    updateSelectizeInput(session, 'clusters', choices = names(data_in_clust()))
  })


  cl_final = reactive({
    req(input$clusters)
    x <- unique(as.vector(data_in_clust_analysis()$Clusters))
    x <- x[!is.na(x)]
    x <- x[x != ""]
    return(x)
  })

  cl_final_analysis = reactive({
    req(input$clusters)
    x <- unique(as.vector(data_in_clust_analysis()$Clusters))
    x <- x[!is.na(x)]
    x <- x[x != ""]
    x <- x[x %in% input$clust_include]
    return(as.vector(x))
  })

  observeEvent (input$confirm, {
    x <- unique(as.vector(data_in_clust()[,input$clusters]))
    x <- x[!is.na(x)]
    x <- x[x != ""]

    # Can use character(0) to remove all choices
    if (is.null(x)) {
      x <- character(0)}

    if (length(x)<=20){
      # Can also set the label and select items
      updateCheckboxGroupInput(session, "clust_include",
                               label = "Which clusters do you want to compare?",
                               choices = x,
                               selected = x
      )}
    else {
      updateCheckboxGroupInput(session, "clust_include",
                               label = "Which clusters do you want to compare?\nWarning: you have more that 20 clusters in your file (max 20)",
                               choices = x[1:20],
                               selected = x[1:20])
    }

  })




  ## keep track of elements inserted and not yet removed
  inserted <- c()

  observeEvent(input$addCluster, {
    btn <- input$addCluster
    id <- paste0('txt', btn)
    insertUI(
      selector = '#placeholder',
      ui = textInput(id, paste("Cluster",btn),
                     value = "")
    )
    inserted <<- c(id, inserted)
  })

  observeEvent(input$reset_cl, {
    ### TODO ###
  })



  ck_DE = eventReactive(input$submit_GO_cl1, {

    #showNotification("Analysis started!")

    #require(org.Hs.eg.db, lib.loc = Rlib)
    #require(org.Mm.eg.db, lib.loc = Rlib)
    #require(org.Dm.eg.db, lib.loc = Rlib)

    ### Define up and down-regulated sets
    df = data_plot()
    up_genes = as.vector(df$GeneID[df$sig == "Significant" & df$log2FC > 0])
    down_genes = as.vector(df$GeneID[df$sig == "Significant" & df$log2FC < 0])

    ### translate across IDs for each cluster
    eg_1 = bitr(up_genes, fromType=input$id_type, toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=input$model_organism)
    eg_2 = bitr(down_genes, fromType=input$id_type, toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=input$model_organism)
    list_genes = list(Up = eg_1$ENTREZID, Down = eg_2$ENTREZID)

    ck <- compareCluster(geneCluster = list_genes,
                         fun = "enrichGO",
                         OrgDb=input$model_organism_cl,
                         ont=input$ont_cl,
                         pAdjustMethod = input$pAdjustMethod,
                         qvalueCutoff = 0.1,
                         pvalueCutoff = 0.1,
                         readable      = TRUE)

      ### return ck object for visualization
      return(ck)


    #showNotification("Analysis is completed!")

  })



  output$dotplot_cl_DE = renderPlot({
    req(input$file)
    dotplot(ck_DE(), showCategory=input$num_show_cl_DE)
  })

  output$output_table_cl_DE = renderDataTable({
    as.data.frame(ck_DE())
  })

  output$downloadData_cl_DE <- downloadHandler(
    filename = function() {
      paste("Results_GO-Analysis_CompareClusters_",Sys.Date(),".tsv", sep="")
    },
    content = function(file) {
      write.csv(as.data.frame(ck_DE()), file, quote = F, sep="\t")
    }
  )


  ck_Tab = eventReactive(input$submit_GO_cl2, {
    
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)
    require(org.Dm.eg.db)

    clust_tab_def = list()
    clust_tab = clusts_EnrichGO()

    for (i in 1:length(clust_tab)){
      v = bitr(clust_tab[[i]], fromType=input$id_type_cl, toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb=input$model_organism_cl)
      clust_tab_def[[i]] = v$ENTREZID
      print(head(clust_tab_def[[i]]))
    }
    names(clust_tab_def) = names(clust_tab)

    ck_Tab_def <- compareCluster(geneCluster = clust_tab_def,
                         fun = "enrichGO",
                         OrgDb=input$model_organism_cl,
                         ont=input$ont_cl,
                         pAdjustMethod = input$pAdjustMethod,
                         qvalueCutoff = 0.1,
                         pvalueCutoff = 0.1,
                         readable      = TRUE)

    return(ck_Tab_def)

  })

  output$comp_cl_Tab = renderPlot(
    dotplot(ck_Tab(), showCategory=input$num_show_cl_Tab)
  )


  output$table_cl_Tab_result = renderDataTable({
    as.data.frame(ck_Tab())
  })


  output$downloadData_cl_Tab <- downloadHandler(
    filename = function() {
      paste("Results_GO-Analysis_CompareClusters_",Sys.Date(),".tsv", sep="")
    },
    content = function(file) {
      write.csv(as.data.frame(ck_Tab()), file, quote = F, sep="\t")
    }
  )


  #####################
  #### sessionInfo ####
  #####################

  output$sessionInfo <- renderText({
    paste(capture.output(sessionInfo()),collapse = "\n")
  })


  #######################
  #### Documentation ####
  #######################

  getPage<-function() {
    return(includeHTML("/data/manke/group/shiny/ferrari/clusterProfiler_GOenrich/ShinyApp_documentation/Documentation_ShinyApp_clusterProfiler.html"))
  }
  
  output$documentation<-renderUI({getPage()})
  #output$documentation<-renderUI({includeHTML("/data/manke/group/shiny/ferrari/clusterProfiler_GOenrich/ShinyApp_documentation/Documentation_ShinyApp_clusterProfiler.html")})
  
  df_sample = reactive({
    df = fread("/data/processing/clusterProfiler_example_data/DEseq_basic_DEresults.tsv")
    return(df)
    })

  
  output$DDS <- downloadHandler(
    filename = function() {
      paste("SampleDataset_",Sys.Date(),".tsv", sep="")
    },
    content = function(file) {
      write.csv(as.data.frame(df_sample()), file, quote = F, row.names = F,sep="\t")
    }
  )
  
  
  output$DDS_v <- downloadHandler(
    filename = "clusterProfiler_vignette.html",
    content = function(con) {
      file.copy(from="/data/manke/group/shiny/ferrari/clusterProfiler_GOenrich/ShinyApp_documentation/Documentation_ShinyApp_clusterProfiler.html", to=con, overwrite =TRUE)
    }
  )
 
  
  
  }

shinyApp(ui, server)
