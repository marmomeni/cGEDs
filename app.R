library(psych)
library(tidyverse)
library(shiny)
library(DT)
library(vroom)
ds<-vroom::vroom("www/Drug sensitivity data (GDSC1).csv")#Drug Sensitivity Data
ex<-vroom::vroom("www/Gene expression data (GDSC).csv ")#Gene Expression 
ui <- fluidPage(
  
  # Application title
  titlePanel("Cancer Drug sensitivity Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("cancer","Select a cancer type",
                  choices=c("Brain lower grade glioma (LGG)","Kidney renal clear cell carcinoma (KIRC)","Esophageal carcinoma (ESCA)",
                  "Breast invasive carcinoma (BRCA)","Stomach adenocarcinoma (STAD)","Mesothelioma (MESO)",
                  "Skin cutaneous melanoma (SKCM)","Lung adenocarcinoma (LUAD)","Glioblastoma multiforme (GBM)",
                  "Head and neck squamous cell carcinoma(HNSC)","Liver hepatocellular carcinoma (LIHC)",
                  "Small cell lung cancer (SCLC)","Neuroblastoma (NB)","Ovarian serous cystadenocarcinoma (OV)",
                  "Colon and rectum adenocarcinoa (COAD/READ) (COREAD)","Multiple myeloma (MM)","Lung squamous cell carcinoma (LUSC)",
                  "Uterine corpus endometrial carcinoma (UCEC)","Pancreatic adenocarcinoma (PAAD)","Acute lymphoblastic leukemia (ALL)",
                  "Head and neck squamous cell carcinoma (HNSC)","Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)",
                  "Medulloblastoma (MB)","Chronic myelogenous leukemia (LCML)","Thyroid carcinoma (THCA)",
                  "Bladder urothelial carcinoma (BLCA)","Prostate adenocarcinoma (PRAD)","Adrenocortical carcinoma (ACC)",
                  " Chronic lymphocytic leukemia (CLL)","Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)",
                  "Acute myeloid leukemia (LAML)",selected=NULL)),
      
      selectizeInput("Genes", "Please enter you desiered genes",choices = colnames(ex[,2:6]),multiple=TRUE),
      
      actionButton("cal","Calculate Correlations")
  ),
    mainPanel(
      
      DT::DTOutput("cortabs"),
      
      downloadButton("download","Download .tsv")
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

  # Provide a vector of drug names
  drugs <- reactive(unique(df()$Drug.name))
  
  # Calculate the correlations and FDRs for each drug separately
  corrs <- NULL
  
  correlations<-eventReactive(input$cal,{
  
    for (i in 1:length(drugs())) {
    
    #Filter the rows related to each drug 
    drug_df <- df() %>%
      filter(Drug.name == drugs()[i])
   
    #drug_df[,4:length(drug_df) refers to gene expression data and drug_df[,3] refers to IC50 column  
    drug_corr <- suppressWarnings(corr.test(drug_df[,4:length(drug_df)], drug_df[,3], method = "pearson",adjust="fdr"))
    
    new_entry <- data.frame(Corr=drug_corr$r, FDR=drug_corr$p.adj) %>%
      mutate(Drug=drugs()[i])
    new_entry$Gene <- row.names(new_entry)
    row.names(new_entry) <- NULL
    corrs <- rbind(corrs, new_entry)
    return(corrs)
  } 
    
  })
  
  # Display the correlation table
  output$cortabs<-DT::renderDT(
    correlations()
    )
  
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


