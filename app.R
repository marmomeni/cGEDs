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

dsGDSC1<-vroom::vroom("www/Drug-sensitivity-data-GDSC1.csv")
dsGDSC2<-vroom::vroom("www/Drug-sensitivity-data-GDSC2.csv")
ex<-vroom::vroom("www/Gene-expression-data-GDSC.csv ")#Gene Expression 


ui <- dashboardPage(
  header <- dashboardHeader(title = dashboardBrand(title = "CGDS app")
  ),
  sidebar <- dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home Page", tabName = "introduction",icon = icon("home")),
      menuItem("Data Selelection & Correlation Calculation", tabName = "dataSelection",icon = icon('mouse-pointer')),
      menuItem("Bubble Plot", tabName = "bubblePlot",icon = icon('chart-line')),
      menuItem("Apply thresholds & Scatter/Boxplot", tabName = "applyThresholds",icon=icon('sliders-h')),
      menuItem("Tutorial", tabName = "tutorial",icon = icon('file-video')),
      menuItem("FAQs", icon = icon("question-circle"), tabName = "faq"),
      menuItem("Contact", tabName = "contact", icon = icon("users")),
      menuItem("Meet Our Team", tabName = "meetteam")
    )
  ),
  body <- dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "introduction",
              fluidRow(
                column(12,
                       align="center",
                       #Introduction block 
                       jumbotron(
                         status = "success",btnName = NULL,
                         title = HTML("<b> CGDS (Cancer Gene-expression Drug-Sensitivity) app </b>"),
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
              fluidRow(
                column(6, 
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
                                          choices = colnames(ex[,3:8]),multiple=TRUE),

                          actionButton("cal","Calculate Correlations", status="success")
                ),
                 column(6,
                          DT::DTOutput("cortabs"),
                          uiOutput("download"), 
                          br(),
                          #uiOutput("threshtab")
                 )
              ),
    ),
    
    tabItem(tabName = "bubblePlot",
            fluidRow(
              column(4,
                     box('This is the Bubble plot page', title = "Bubble plot",  
                         status = "primary", solidHeader = TRUE,
                         collapsed = FALSE, width=12)                    
              )
            )
    ),

    tabItem(tabName = "applyThresholds",
        fluidPage(
              column(12,
                     div(style = "display:inline-block; float:left",
                         actionButton('to_introduction', label = 'Home', status = "success")),
                     div(style = "display:inline-block; float:right",
                         actionButton('backto_tutorial', label = 'See Tutorial', status = "success"))
              ),
              column(12, align="center",
                     HTML("<h5>Apply thresholds to select a Gene/drug pair for Scatter/Boxplot</h5>")
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
                   actionBttn("Thre","Apply Thresholds",style="pill",color="success",size = "sm")))
            ,column(5,align="center",offset = 1,
                    wellPanel(DT::DTOutput("Sigcors")), 
             )
          ),
          hr(),
          fluidRow(column(6,align="center",
                   br(),
                   uiOutput("selGenedrug"),
                   br(),
                   br(),
                   prettyCheckbox("scatterLabel","Show cell line names",value = TRUE
                                  ,status = "success", outline = TRUE),
                   prettyCheckbox("ShowBoxplot","Show marginal boxplots",value = TRUE,
                                  status = "success", outline = TRUE),
                   div(id = "Col0"),
                   div(id = "Col1"),
                   div(id = "Col2")
            ),
            column(4,align="center",
                   br(),
                   br(),
                   plotOutput("scatterplt",width = "100%")
            )
          )
  ),
  tabItem(tabName = "bubblePlot",
          fluidRow(
            column(4,
                   box('This is the Bubble plot page', title = "Bubble plot",  
                       status = "primary", solidHeader = TRUE,
                       collapsed = FALSE, width=12)                    
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

  observeEvent(input$to_dataSelect, {
    updateTabItems(session, "tabs", selected = "dataSelection")
  }
  )
  observeEvent(input$to_tutorial, {
    updateTabItems(session, "tabs", selected ="tutorial")
  }
  )
  
  observeEvent(input$backto_tutorial, {
    updateTabItems(session,"tabs", selected ="tutorial")
  }
  )
  observeEvent(input$to_faq, {
    updateTabItems(session, "tabs", selected ="faq")
  }
  )
  
  observeEvent(input$to_introduction, {
    updateTabItems(session, "tabs",selected = "introduction")
  })
  
  # DATA SELECT BUTTONS
  observeEvent(input$backTo_introduction, {
    updateTabItems(session, "tabs", "introduction")
  }
  )
  observeEvent(input$to_demo, {
    updateTabItems(session, "tabs", "demo")
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
    drugs <- unique(df()$Drug.name)
    corrs <- NULL
    
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
  
  # "Apply thresholds and visualization" button appears when clicking on the "Calculate Correlations" button 
  observeEvent(input$cal, {
    output$vistab <- renderUI({actionBttn("vistab","Next"
                                          ,style="pill",color="success",size = "sm")})
  })
  
  #When clicking on the "Next" button, the tab switches to
  #the apply threshold tab
  observeEvent(input$vistab, {
    updateTabItems(session, "tabs", selected = "applyThresholds")
  })
  
 #Apply thresholds
  sigcors<-eventReactive(input$Thre,{
  sigcors1<-subset(correlations(), FDR< input$FDRThr & Corr> min(input$PosCorThre))
  sigcors2<-subset(correlations(), FDR< input$FDRThr & Corr< max(input$NegCorThre))
  sigcors<-(rbind(sigcors1,sigcors2))
  })
  
  #Provide Gene/Drug pair list for the Drop-down of choosing Gene/Drug pair for the visualization
  sigcors4<-reactive({
  GeneDrug <-paste(sigcors()$Gene ," / ", sigcors()$Drug)
  sigcors3<-cbind(sigcors(),GeneDrug)
  })
  
  #Drop-down for choosing Gene/Drug pair for the visualization
   observeEvent(input$Thre,{output$selGenedrug<-renderUI({
    selectInput("selGenedrug", "Please choose desiered Gene/Drug pair for the visualization",
                   choices = sigcors4()$GeneDrug,multiple=FALSE)
  }) })

  #Scatter/boxplot
   Scatter<-reactive({
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

   # Display the correlation table
  output$cortabs<-DT::renderDT(
    correlations()
  )
  output$Sigcors<-DT::renderDT(
    sigcors()
  )
  #output$volcanopl<-renderPlot(
  #  Volcano()
  #)
  output$scatterplt<-renderPlot(
    Scatter(),res = 96, height = 600, width = 600 
   )
             
}

shinyApp(ui = ui, server = server)


