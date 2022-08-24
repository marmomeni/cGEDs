library(shinydashboard)
library(shiny)
library(bs4Dash)
library(fresh)
library(psych)
library(tidyverse)
library(DT)
library(vroom)
ds<-vroom::vroom("www/Drug sensitivity data (GDSC1).csv")#Drug Sensitivity Data
ex<-vroom::vroom("www/Gene expression data (GDSC).csv ")#Gene Expression 


ui <- dashboardPage(
  header <- dashboardHeader(title = dashboardBrand(title = "CGDS app")
  ),
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Home Page", tabName = "introduction"),
      menuItem("Select Data", tabName = "dataSelect"),
      menuItem("See Demo", tabName = "demo"),
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
                           actionButton('to_demo', label = 'See Demo', status = "success"),
                           actionButton('to_faq', label = 'FAQ', status = "success"))
                ),
              )
      ),
    tabItem(tabName = "dataSelect",fluidPage(
       navbarPage("CGDS (Cancer Gene-expression Drug-Sensitivity app)",id="inTabset",
          tabPanel("Correlation Calculation",
            sidebarPanel(
              selectInput("cancer","Select a cancer type",
                  choices=c("Brain lower grade glioma (LGG)",
                  "Kidney renal clear cell carcinoma (KIRC)",
                  "Esophageal carcinoma (ESCA)",
                  "Breast invasive carcinoma (BRCA)" ,
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
                  "Acute myeloid leukemia (LAML)",selected=NULL)),
              selectizeInput("Genes", "Please enter your desiered genes",
                     choices = colnames(ex[,3:10]),multiple=TRUE),
              actionButton("cal","Calculate Correlations",class="btn btn-success"),
              br(),
           ),
 
            mainPanel(
      
              DT::DTOutput("cortabs"),
              uiOutput("download"),
              br(),
              uiOutput("vistab")
              
            )
           ),
      
          #selectInput("gene-drug","Please select the gene-drug assossiation you want to visualize",choices=)
          tabPanel("Visualization",
            sidebarPanel(
              numericInput("FDRThr","Choose Gene/drug pairs with FDRs less than:", value = 0.05),
              sliderInput("PosCorThre", "Choose Gene/drug pairs with correlations more than:",min = 0, max =1,value = 0.7,step = 0.1),
              sliderInput("NegCorThre", "Choose Gene/drug pairs with correlations less than:",min = -1, max =0,value = -0.7,step = 0.1),
              br(),
              actionButton("Thre","Apply Thresholds",class="btn btn-success")
            ),
              br(),
              br(),
            mainPanel(
              DT::DTOutput("Sigcors"),
              br(),
              br(),
              uiOutput("outputUI"),
              plotOutput("scatterplt")
            )
           )
      )
  )
    ),
  tabItem(tabName = "demo",
          fluidRow(
            column(4,
                   box('This is the demo page', title = "Demo Page",  
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
  
  # Cancer type selection by the user
  ds2 <- reactive(ds %>% 
    filter(ds$`Cancer-Type`== input$cancer))

  # remove the cancer type column since it's not necessary anymore
  ds3<-reactive(ds2()[-4])
  
  # Gene selection by the user
  ex2<-reactive(ex[,input$Genes])
  
  # Add Cell line column to ex2 needed for merging step
  ex3 <- reactive(cbind('Cell line' =ex[1], ex2()))
  
  # Merge the two tables
  df <- reactive(merge(x = ds3(), y = ex3(), by ="Cell line"))
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
    output$vistab <- renderUI({actionButton("vistab" ,"Apply thresholds and visualization",class="btn btn-success")})
  })
  
  #When clicking on the "Apply thresholds and visualization button appears" button, the tab switches to
  #the Visualization tab
  observeEvent(input$vistab, {
    updateTabsetPanel(session, "inTabset",selected = "Visualization")
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
   observeEvent(input$Thre,{output$outputUI<-renderUI({
    selectInput("outputUI", "Please choose desiered Gene/Drug pair for the visualization",
                   choices = sigcors4()$GeneDrug,multiple=FALSE)
  }) })
  
   #Volcano Plot
   
   #Volcano<-reactive({
   #drug_example <- subset(correlations(), Drug == sigcors4()[which(sigcors4()$GeneDrug == input$outputUI) , 3])
   
   #p <- ggplot(drug_example, aes(x = drug_example$Corr , y = -log10(drug_example$FDR))) + geom_point() + theme_minimal() 
   
   #p + geom_vline(xintercept=c(0), col="red") #+ ggtitle(past(correlations()[which(correlations()$GeneDrug == input$outputUI) , 3],"drug"))
   
   #})
  
   Scatter<-reactive({
     drug_example <- subset(df(), Drug.name == sigcors4()[which(sigcors4()$GeneDrug ==input$outputUI) , 3])
    
     gene<-as.character(sigcors4()[which(sigcors4()$GeneDrug ==input$outputUI), 4])
     ggplot(drug_example,aes(drug_example[,gene],
                             IC50,
                             label= `Cell line`))+geom_point(size=2)+
                             geom_text(nudge_x = 0, nudge_y = 0.2,size=6,color="darkcyan")+
                             theme_bw()+theme(text = element_text(size=20))+geom_smooth(method=("lm"))
    
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
    Scatter(),res = 96 
   )
             
}

shinyApp(ui = ui, server = server)


