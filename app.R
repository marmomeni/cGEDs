library(psych)
library(tidyverse)
library(shiny)
library(DT)
library(vroom)
ds<-vroom::vroom("www/Drug sensitivity data (GDSC1).csv")#Drug Sensitivity Data
ex<-vroom::vroom("www/Gene expression data (GDSC).csv ")#Gene Expression 
ui <- fluidPage(
  
  # Application title
  navbarPage(
    "CGDS (Cancer Gene-expression Drug Sensitivity app)",
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
      br(),
      br(),
      selectizeInput("Genes", "Please enter your desiered genes",
                     choices = colnames(ex[,3:8]),multiple=TRUE),
      br(),
      br(),
      actionButton("cal","Calculate Correlations")
      ),
 
    mainPanel(
      
      DT::DTOutput("cortabs"),
      
      downloadButton("download","Download .tsv"),

    )
  ),
      
      #selectInput("gene-drug","Please select the gene-drug assossiation you want to visualize",choices=)
  tabPanel("Visualization",
    sidebarPanel(
      numericInput("FDRThr","Choose Gene/drug pairs with FDRs less than:", value = 0.05),
      br(),
      br(),
      sliderInput("PosCorThre", "Choose Gene/drug pairs with correlations more than:",min = 0, max =1,value = 0.7,step = 0.1),
      sliderInput("NegCorThre", "Choose Gene/drug pairs with correlations less than:",min = -1, max =0,value = -0.7,step = 0.1),
      br(),
      br(),
      actionButton("Thre","Apply Thresholds")
      ),
    mainPanel(
      DT::DTOutput("Sigcors"),
      br(),
      br(),
      uiOutput("outputUI"),
      #plotOutput("volcanopl")
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
  
  # Add the ability to download the correlation table
  output$download <- downloadHandler(
    filename = function() {
      "Pearson Correlations and FDRs.tsv"
    },
    content = function(file) {
      vroom::vroom_write(correlations(), file)
    }
  )                 
}

shinyApp(ui = ui, server = server)


