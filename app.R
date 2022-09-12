library(shinydashboard)
library(shiny)
library(bs4Dash)
library(fresh)
library(psych)
library(tidyverse)
library(DT)
library(vroom)
library(ggrepel)
library(ggExtra)
library(colourpicker)
library(shinyWidgets)
library(shinyjs)
library(dplyr)


dsGDSC1<-vroom::vroom("www/Drug-sensitivity-data-GDSC1.csv")
dsGDSC2<-vroom::vroom("www/Drug-sensitivity-data-GDSC2.csv")
ex<-vroom::vroom("www/Gene-expression-data-GDSC.csv ")#Gene Expression 


ui <- dashboardPage(
  header <- dashboardHeader(title = dashboardBrand(title = "cGEDs app")
  ),
  sidebar <- dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home Page", tabName = "introduction",icon = icon("home")),
      menuItem(text=tags$div("Data Selelection &",tags$br(), "Correlation Calculation",style= "display: inline-block;vertical-align:middle"), tabName = "dataSelection",icon = icon('mouse-pointer')),
      menuItem("Apply thresholds", tabName = "applyThresholds",icon=icon('sliders-h')),
      menuItem("Bubble Plot", tabName = "bubblePlot",icon = icon('chart-line')),
      menuItem("Scatter/Boxplot", tabName = "scatterBoxplot",icon = icon('chart-line')),
      menuItem("Tutorial", tabName = "tutorial",icon = icon('file-video')),
      menuItem("FAQs", icon = icon("question-circle"), tabName = "faq"),
      menuItem("Contact", tabName = "contact", icon = icon("users")),
      menuItem("Meet Our Team", tabName = "meetteam")
    )
  ),
  body <- dashboardBody(
    tabItems(
      # Home page content
      tabItem(tabName = "introduction",
              fluidRow(
                column(12,
                       align="center",
                       #Introduction block 
                       jumbotron(
                         status = "success",btnName = NULL,
                         title = HTML("<b> cGEDs (cancer Gene-Expression Drug-sensitivity) app </b>"),
                         lead = "An application for finding drug effectivity biomarkers in different cancer types",
                         "This app is created by STEM-Away RShiny Project Team - Session 1, 2022",
                         href = "https://stemaway.com/" 
                       ),
                       div(style ="display:inline-block", 
                           actionButton('to_dataSelect', label = 'Begin', status = "success"),
                           actionButton('to_tutorial', label = 'Tutorial', status = "success"),
                           actionButton('to_faq', label = 'FAQ', status = "success")
                       )
                 )
              ),
      ),

      tabItem(tabName = "dataSelection",
         fluidPage(
                column(12,
                       div(style = "display:inline-block; float:left",
                           actionButton('to_introduction', label = 'Home', status = "success")),
                       div(style = "display:inline-block; float:right",
                           actionButton('to_tutorial', label = 'See Tutorial', status = "success"))
                ),
                column(12, align="center",
                       HTML("<h5>Choose from the options given to begin</h5>")
                ),
                hr(),
          ),          
          fluidRow(
                column(4,align="center",offset = 1,
                         selectInput("dataset","Select a drug sensitivity and gene expression dataset",
                                choices=c("GDSC1","GDSC2"),selected=NULL),
                         selectInput("cancer","Select a cancer type",
                                choices=c("Brain lower grade glioma (LGG)",
                                          "Kidney renal clear cell carcinoma (KIRC)",
                                          "Esophageal carcinoma (ESCA)",
                                          "Breast invasive carcinoma (BRCA)",
                                          "Stomach adenocarcinoma (STAD)","Mesothelioma (MESO)",
                                          "Skin cutaneous melanoma (SKCM)","Lung adenocarcinoma (LUAD)",
                                          "Glioblastoma multiforme (GBM)",
                                          "Head and neck squamous cell carcinoma(HNSC)",
                                          "Liver hepatocellular carcinoma (LIHC)",
                                          "Small cell lung cancer (SCLC)","Neuroblastoma (NB)",
                                          "Ovarian serous cystadenocarcinoma (OV)",
                                          "Colon and rectum adenocarcinoa (COAD/READ) (COREAD)",
                                          "Multiple myeloma (MM)",
                                          "Lung squamous cell carcinoma (LUSC)",
                                          "Uterine corpus endometrial carcinoma (UCEC)",
                                          "Pancreatic adenocarcinoma (PAAD)",
                                          "Acute lymphoblastic leukemia (ALL)",
                                          "Head and neck squamous cell carcinoma (HNSC)",
                                          "Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)",
                                          "Medulloblastoma (MB)","Chronic myelogenous leukemia (LCML)",
                                          "Thyroid carcinoma (THCA)",
                                          "Bladder urothelial carcinoma (BLCA)","Prostate adenocarcinoma (PRAD)",
                                          "Adrenocortical carcinoma (ACC)"," Chronic lymphocytic leukemia (CLL)",
                                          "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)",
                                          "Acute myeloid leukemia (LAML)"),selected=NULL)
                          ,
                          selectizeInput("Genes", "Please enter your desiered genes",
                                          choices = colnames(ex[,3:50]),multiple=TRUE),
                          useSweetAlert(),
                          actionButton("cal","Calculate Correlations", status="success")
                ),
                 column(5,align="center",offset = 1,wellPanel(
                          DT::DTOutput("cortabs"),
                          br(),
                          uiOutput("download")
                          #uiOutput("threshtab")
                          )
                 )
              ),
          fluidRow(
                column(12,
                         hr()
                )
          ),
         fluidRow(
           column(12,align="center",
                  uiOutput("to_next")
           )
         )         
    ),
    tabItem(tabName = "applyThresholds",
            fluidPage(
              column(12,
                     div(style = "display:inline-block; float:left",
                         actionButton('backto_introduction', label = 'Home', status = "success")),
                     div(style = "display:inline-block; float:right",
                         actionButton('backto_tutorial', label = 'See Tutorial', status = "success"))
              ),
              column(12, align="center",
                     HTML("<h5>Apply thresholds to select the most associated Gene/drug pairs for the visualization</h5>")
              ),
              hr(),
            ),      
            fluidRow(
              column(4,align="center",offset = 1,
                     numericInput("FDRThr","Choose Gene/drug pairs with FDRs less than:", value = 0.05),
                     wellPanel(
                       sliderInput("PosCorThre", "Choose Gene/drug pairs with correlations more than:",min = 0, max =1,value = 0.7,step = 0.1),
                       sliderInput("NegCorThre", "Choose Gene/drug pairs with correlations less than:",min = -1, max =0,value = -0.7,step = 0.1),
                       br(),
                       actionBttn("Thre","Apply Thresholds",style="pill",color="success",size = "sm"))),
              column(5,align="center",offset = 1,
                      wellPanel(DT::DTOutput("Sigcors")), 
              )
            ),
            fluidRow(
              column(12,
                     hr()
              )
            ),
            fluidPage(
              column(12,align="center",
                     br(),
                     div(style ="display:inline-block", 
                         uiOutput('to_scatterPlot')),
                     div(style ="display:inline-block", 
                         uiOutput('to_bubblePlot'))
              )
              
            )
    ), 
            
    tabItem(tabName = "bubblePlot",
            fluidRow(
              column(12, align="center",
                     HTML("<h5>Bubble Plot,  </h5>")
              ),hr(),
            ),
              
            fluidRow(
              column(12,align="center",
                     actionButton('plotBubbleplot', label = 'Plot Bubble Plot', status = "success"),
                     plotOutput("bubble",width = "auto",height = "auto"),
                     hr(),
                     uiOutput("bubbledownload", label = "Download")
              )
            )
    ),

    tabItem(tabName = "scatterBoxplot",
          fluidRow(
            column(12, align="center",
                   HTML("<h5>Select among associated gene/drug pairs to be visualized by a scatter plot with a marginal boxplot </h5>")
            )
          ),
            hr(),
          fluidRow(
            column(6,align="center",
                   br(),
                   uiOutput("selGenedrug"),
                   br(),
                   br(),
                   uiOutput("scatterLabel"),
                   uiOutput("ShowBoxplot"),
                   div(id = "Col0"),
                   div(id = "Col1"),
                   div(id = "Col2")
            ),
            column(4,align="center",
                   br(),
                   br(),
                   plotOutput("scatterplt",width = "100%"),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   uiOutput("scatterdownload", label = "Download")
            )
          ),
          fluidRow(column(12,align="center",
                   #uiOutput('downloadScatter')
                   hr()
                   )      
          )
          
    ),

  tabItem(tabName = "tutorial",
          fluidRow(
            column(4,
                   box('This is the tutorial page', title = "tutorial",  
                       status = "primary", solidHeader = TRUE,
                       collapsed = FALSE, width=12)                    
            )
          )
  ),
  tabItem(tabName = "faq",
          fluidRow(
            column(4,
                   box('Questions here', title = "FAQ Page",  
                       status = "primary", solidHeader = TRUE, collapsible = T,
                       collapsed = FALSE, width=12)                    
            )
          )
  ),
  tabItem(tabName = "contact",
          fluidRow(
            column(4,
                   box('contact', title = "Contact us",  
                       status = "primary", solidHeader = TRUE, collapsible = T,
                       collapsed = FALSE, width=12)                    
            )
          )
  ),
  tabItem(tabName = "meetteam",
          fluidRow(
            column(4,
                   box('content goes here', title = "The team",  
                       status = "primary", solidHeader = TRUE, collapsible = T,
                       collapsed = FALSE, width=12)                    
            )
          )
       )
    )
  )
)
server <- function(input, output,session) {
  
  # Home page buttons
  observeEvent(input$to_dataSelect, {
    updateTabItems(session, "tabs", selected = "dataSelection")
  }
  )
  observeEvent(input$to_tutorial, {
    updateTabItems(session, "tabs", selected ="tutorial")
  }
  )
  
  observeEvent(input$to_faq, {
    updateTabItems(session, "tabs", selected ="faq")
  }
  )
  
  # Data Selection & Correlation Calculation Buttons
  observeEvent(input$to_introduction, {
    updateTabItems(session, "tabs", "dataSelection")
  }
  )
  observeEvent(input$to_tutorial, {
    updateTabItems(session, "tabs", "tutorial")
  }
  )
  observeEvent(input$to_next, {
    updateTabItems(session, "tabs", "applyThresholds")
  }
  )  

  
  #Apply Thresholds button
  observeEvent(input$to_bubblePlot, {
    updateTabItems(session, "tabs", "bubblePlot")
  }
  )
  
  observeEvent(input$to_scatterPlot, {
    updateTabItems(session, "tabs", "scatterBoxplot")
  }
  )  
  # Bubble Plot page buttons
  
  # Apply thresholds & Scatter/boxplot buttons
  observeEvent(input$to_introduction, {
    updateTabItems(session, "tabs", "backto_introduction")
  }
  )
  observeEvent(input$to_tutorial, {
    updateTabItems(session, "tabs", "backto_tutorial")
  }
  )   
  
# Dataset and cancer type selection by the user
  dataselect<-reactive({
  if (input$dataset=="GDSC1"){
  ds <- dsGDSC1 %>% 
    filter(dsGDSC1$`Cancer-Type`== input$cancer)
  }
  else if(input$dataset=="GDSC2"){
  ds<-dsGDSC2 %>% 
    filter(dsGDSC2$`Cancer-Type`== input$cancer)
  } 
  
  })
  # remove the cancer type column since it's not necessary anymore
  ds2<-reactive(dataselect()[-4])
  
  # Gene selection by the user
  ex2<-reactive(ex[,input$Genes])
  
  # Add Cell line column to ex2 needed for merging step
  ex3 <- reactive(cbind('Cell line' =ex[1], ex2()))
  
  # Merge the two tables
  df <- reactive(merge(x = ds2(), y = ex3(), by ="Cell line"))
  correlations<-eventReactive(input$cal,{
    
    # Provide a vector of drug names
    drugs <- df() %>%
      select(Drug.name, 'Cell line') %>%
      group_by(Drug.name) %>%
      summarise(Num_cell_lines=n()) %>%
      filter(Num_cell_lines > 2)
    
    drugs <- unique(df()$Drug.name)
    corrs <- NULL
    
    # Progress bar code
    progressSweetAlert(
      session = session, id = "myprogress",
      title = "Work in progress",
      display_pct = TRUE, value = 0
    )
    # Calculate the correlations and FDRs for each drug separately
    
    for (i in 1:length(drugs)) {
    
    #Filter the rows related to each drug 
      drug_df <- df() %>%
        filter(Drug.name == drugs[i])
   
    #drug_df[,5:ncol(drug_df) refers to gene expression data and drug_df[,3]
    #refers to IC50 column  
      drug_corr <- suppressWarnings(corr.test(drug_df[,5:length(drug_df)], drug_df[,3]
      , method = "pearson",adjust="fdr"))
    
      new_entry <- data.frame(Corr=drug_corr$r, FDR=drug_corr$p.adj) %>%
        mutate(Drug=drugs[i])
      new_entry$Gene <- row.names(new_entry)
      row.names(new_entry) <- NULL
      corrs <- rbind(corrs, new_entry)
    }
    for (i in seq_len(50)) {
      Sys.sleep(0.1)
      updateProgressBar(
        session = session,
        id = "myprogress",
        value = i*2)
    }
    closeSweetAlert(session = session)
    sendSweetAlert(
      session = session,
      title =" Correlation Calculation completed !",
      type = "success"
    )
    return(corrs)
  })
  # Add the ability to download the correlation table
  # Download button appears after clicking on the calculate button using observeEvent and renderUI
  observeEvent(input$cal, {
    output$download <- renderUI({ 
      downloadHandler(
        filename = function() {
          "Pearson Correlations and FDRs.tsv"
        },
        content = function(file) {
          vroom::vroom_write(correlations(), file)
        }
        
      ) 
    }) 
  })
 
  # "Next" button appears when clicking on the "Correlation Calculation" button
  observeEvent(input$cal, {output$to_next<-renderUI({actionButton("to_next","Next",status="success")})
  })
  
  # "Bubble Plot" and "Scatter/Boxplot" buttons appear when clicking on the "Apply Thresholds" button 
  observeEvent(input$Thre, {output$to_bubblePlot<-renderUI({actionButton("to_bubblePlot","Bubble Plot",status="success")})
  })
  
  observeEvent(input$Thre, {output$to_scatterPlot<-renderUI({actionButton("to_scatterPlot","Scatter/Boxplot",status="success")})
  })
  
  

 #Apply thresholds
  sigcors<-eventReactive(input$Thre,{
    sigcors1<-subset(correlations(), FDR< input$FDRThr & Corr> min(input$PosCorThre))
    sigcors2<-subset(correlations(), FDR< input$FDRThr & Corr< max(input$NegCorThre))
    sigcors<-(rbind(sigcors1,sigcors2))
  })
  
 # Bubble plot
  Bubbleplot<-reactive({
    data <- correlations() %>%
      rename(STARS=Corr, FDR_NUM=FDR) %>%
      separate(STARS, into=c('CORR'), sep='\\*', remove=FALSE, extra='drop') %>%
      mutate(CORR = as.numeric(CORR))
    
    
    bubble_df <- data %>%
      filter(Drug %in% sigcors()$Drug) %>%
      filter(Gene %in% sigcors()$Gene) %>%
      select(Drug, Gene, CORR, FDR_NUM) %>%
      mutate(FDR = ifelse(FDR_NUM <= 0.0001, '<=0.0001', 
                          ifelse(FDR_NUM <= 0.001, '<0.001', 
                                 ifelse(FDR_NUM <= 0.01, '<0.01', 
                                        ifelse(FDR_NUM <= 0.05, '<0.05', '>0.05'))))) %>%
      mutate(FDR = factor(FDR, levels=c('>0.05', '<0.05', '<0.01', '<0.001', '<=0.0001'))) %>%
      mutate(Significant = ifelse(FDR_NUM <= 0.05, '<=0.05', '>0.05')) %>%
      mutate(Significant = factor(Significant, levels=c('>0.05', '<=0.05'))) %>%
      rename(Correlation=CORR)
    
    limit <- max(abs(bubble_df$Correlation)) * c(-1, 1)
    ggplot(bubble_df, aes(x=Drug, y=Gene, alpha=Significant, size=FDR, colour=Correlation))+
      geom_point()+
      scale_colour_distiller(type = "div", palette='RdBu', limit=limit)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
  })
  
 
  
  #Provide Gene/Drug pair list for the Drop-down of choosing Gene/Drug pair for the visualization
  sigcors4<-reactive({
  GeneDrug <-paste(sigcors()$Gene ," / ", sigcors()$Drug)
  sigcors3<-cbind(sigcors(),GeneDrug)
  })
  
  #Drop-down for choosing Gene/Drug pair for the visualization
   observeEvent(input$Thre,{output$selGenedrug<-renderUI({
    selectInput("selGenedrug", "Please choose desiered Gene/Drug pair for the visualization",
                   choices = c("",sigcors4()$GeneDrug),multiple=FALSE,selected=NULL)
  }) })

   observeEvent(req(outpu$scatterplt),{output$scatterLabel<-renderUI({
     prettyCheckbox("scatterLabel","Show cell line names",value = TRUE
                    ,status = "success", outline = TRUE)
   }) })
   
   observeEvent(req(outpu$scatterplt),{output$ShowBoxplot<-renderUI({
     prettyCheckbox("ShowBoxplot","Show marginal boxplots",value = TRUE,
                    status = "success", outline = TRUE)
   }) })

  #Scatter/boxplot
   
   Scatter<-reactive({
     req(input$selGenedrug)
     drug<-sigcors4()[which(sigcors4()$GeneDrug ==input$selGenedrug) , 3]
     drug_df <- subset(df(), Drug.name == drug)
     gene<-as.character(sigcors4()[which(sigcors4()$GeneDrug ==input$selGenedrug), 4])
     med=median(drug_df[,gene])
     drug_df$GeneExpressLevel = ifelse (drug_df[,gene] >= med, "high", "low")
    
     if(input$ShowBoxplot==TRUE){
        if(input$scatterLabel==TRUE){
          x<-ggplot(drug_df,aes(drug_df[,gene],IC50,label=`Cell line`))+
          theme_bw()+theme(text = element_text(size=12), legend.position='bottom')+
          geom_smooth(method=("lm"))+
          labs( x = paste("Expression levels of",gene),y= paste("IC50 of",drug))+
          geom_point(size=2, aes(colour=GeneExpressLevel))+
          scale_colour_manual(values=c(input$col1, input$col2))+
          geom_text_repel()
          ggMarginal(x,type="boxplot",groupColour=TRUE,groupFill = TRUE)
        }
        else{
          x<-ggplot(drug_df,aes(drug_df[,gene],IC50,label=`Cell line`))+
          theme_bw()+theme(text = element_text(size=12), legend.position='bottom')+
          geom_smooth(method=("lm"))+
          labs( x = paste("Expression levels of",gene),y= paste("IC50 of",drug))+
          geom_point(size=2, aes(colour=GeneExpressLevel))+
          scale_colour_manual(values=c(input$col1,input$col2))
          ggMarginal(x,type="boxplot",groupColour=TRUE,groupFill = TRUE)
        }
     }
     else{
        if(input$scatterLabel==FALSE){
          ggplot(drug_df,aes(drug_df[,gene],IC50,label=`Cell line`))+
          theme_bw()+theme(text = element_text(size=12), legend.position='bottom')+
          geom_smooth(method=("lm"))+
          labs( x = paste("Expression levels of",gene),y= paste("IC50 of",drug))+
          geom_point(size=1)
         }
        else{
          ggplot(drug_df,aes(drug_df[,gene],IC50,label=`Cell line`))+
          theme_bw()+theme(text = element_text(size=12), legend.position='bottom')+
          geom_smooth(method=("lm"))+
          labs( x = paste("Expression levels of",gene),y= paste("IC50 of",drug))+
          geom_point(size=1)+geom_text_repel()
         }
     }
   })
   
   # ObserveEvent functions related to selecting the color of scatter/boxplots
   observeEvent(input$ShowBoxplot,{
     if (input$ShowBoxplot==TRUE){
       insertUI(
         selector = "#Col1",
         ui=colourpicker::colourInput("col1", "Select colour of boxplots related to the high expression cell lines"
       ,showColour="both", value = "5158AD")
       )
     }
     else if(input$ShowBoxplot==FALSE){
       removeUI(
         selector = "div#Col1 > div"
       )
     } 
   })
  
   observeEvent(input$ShowBoxplot,{
     if (input$ShowBoxplot==TRUE){
       insertUI(
         selector = "#Col2",
         ui=colourpicker::colourInput("col2", "Select colour of boxplots related to the low expression cell lines"
         ,showColour="both", value = "F78A25")
       )
     }
     else if(input$ShowBoxplot==FALSE){
       removeUI(
         selector = "div#Col2 > div"
       )
     } 
   })
  
   observeEvent(input$ShowBoxplot,{
     if (input$ShowBoxplot==FALSE){
       insertUI(
         selector = "#Col0",
         ui=colourpicker::colourInput("col0", "Select colour of dots"
         ,showColour="both", value = "F78A25")
       )
     }
     else if(input$ShowBoxplot==TRUE){
       removeUI(
         selector = "div#Col0 > div"
       )
     } 
   })
   
   
   output$bubble <-renderPlot(
    Bubbleplot(),res = 96, height =function(){length(unique(sigcors()$Gene))*15+450} , width = function(){length(unique(sigcors()$Drug))*35+300}
  )
   
   
   observeEvent(input$plotBubbleplot, {
     # The download button of the scatter plot 
     output$bubbledownload <-renderUI({ downloadHandler(
       filename =  function() {
         "Bubble Plot.png"
       },
       # content is a function with argument file. content writes the plot to the device
       content = function(file) {
         device <- function(..., width, height) {
           grDevices::png(..., width = width, height = height,
                          res = 300, units = "in")}
         ggsave(file, plot = Bubbleplot(), device = device)
       } 
     )
     })
   })
   
   output$scatterplt<-renderPlot(Scatter(),res = 96, height = 600, width = 600)
   
   # The download button of the scatter plot spears after selGenedrug drop-down works
   observeEvent(input$selGenedrug, {
     # The download button of the scatter plot 
     output$scatterdownload <-renderUI({ downloadHandler(
     filename =  function() {
       "Scatter Plot.png"
     },
     # content is a function with argument file. content writes the plot to the device
     content = function(file) {
       device <- function(..., width, height) {
         grDevices::png(..., width = width, height = height,
                        res = 300, units = "in")}
       ggsave(file, plot = Scatter(), device = device)
     } 
   )
   })
   })
   
   output$cortabs<-DT::renderDT(
     correlations()
  )
  output$Sigcors<-DT::renderDT(
    sigcors()
  )
  
  
  #output$scatterplt<-renderPlot(
  #  Scatter(),res = 96, height = 600, width = 600,
    #ggsave("plot.pdf",Scatter()) 
  # )
 

}

shinyApp(ui = ui, server = server)


