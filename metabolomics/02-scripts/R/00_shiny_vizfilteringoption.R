library(tidyverse)
library(ggplot2)
library(plotly)
library(ggExtra)
library(vegan)
library(htmlwidgets)
library(Hmisc)
library(shiny)

#define UI for app
ui <- fluidPage(
  
  #App title 
  titlePanel("Metabolomics Viewer App"),
  
  #Sidebar layout with input and output definition
  sidebarLayout(
    
    #Sidebar panel for inputs
    sidebarPanel(
      
      selectInput('data_choice', "Select Dataframe to Check",
                  choices = dflist,
                  selected = dflist[1]),
      selectInput("comparison", "Select Linear Model Comparisons", 
                  choices = unique(get(paste0(strsplit(dflist[1], "_")[[1]][1], "_", strsplit(dflist[1], "_")[[1]][3],
                                              "_maaslin_onlyinfl"))$metadata),
                  selected = unique(get(paste0(strsplit(dflist[1], "_")[[1]][1], "_", strsplit(dflist[1], "_")[[1]][3],
                                               "_maaslin_onlyinfl"))$metadata))
      
    ),
    #Main panel for outputs
    mainPanel(
      plotlyOutput("volcano_plot"),
      plotOutput(outputId = "feature_violinplot"),
      uiOutput("feature comparisons")
    )
  )
)

server <- function(input, output){
  
  val_df <- eventReactive(input$data_choice, {get(input$data_choice)})
  
  linmod_df <- eventReactive(input$data_choice, {
    name <- paste0(strsplit(input$data_choice, "_")[[1]][1], "_", strsplit(input$data_choice, "_")[[1]][3],
                   "_maaslin_onlyinfl") #change file of linear model to use here
    get(name)
  })
  
  
  d.filtered <- reactive({linmod_df() %>%
      filter(metadata == input$comparison) %>%
      select(one_of(c("coef", "qval", "feature", "nameFix", "metadata", "value"))) %>%
      mutate(log10_pval_adj = -log10(qval),
             upordown = ifelse(coef > 0, "Up", "Down"),
             log2_FC = coef / log10(2)) %>%
      mutate(upordown = ifelse(log10_pval_adj <1.3, "Not Significant", upordown))
  })
  
  # d.boxplot <- reactive({
  #   val_df
  #   
  #   
  # })
  
  gp <- reactive({ggplotly(
    ggplot(d.filtered(), aes(x = log2_FC, y = log10_pval_adj, color = upordown, text = nameFix)) +
      geom_point() +
      theme_bw() +
      labs(x = "log2FC", y = "-log10(p-value)", color = NULL) + 
      scale_color_manual(values = c("red", "blue", "grey"), 
                         limits = c("Up", "Down", "Not Significant"))
  )})
  
  output$volcano_plot <- renderPlotly(gp())
  
  s <- reactive(event_data("plotly_click"))
  
  feat_boxplot <- eventReactive(s(),{
    feat_nm <- d.filtered() %>%
      filter(log2_FC == s()$x & log10_pval_adj == s()$y) %>%
      select(nameFix) %>% as.character() #select the unique name
    
    m <- d.filtered() %>%
      filter(log2_FC == s()$x & log10_pval_adj == s()$y) %>%
      select(metadata) %>% as.character() #select the unique metadata
    
    data <- val_df() %>% filter(nameFix == feat_nm) %>%
      select(one_of(append(unique(linmod_df()$metadata), c("nameFix", "log_pkh_imputed", "donor_type")))) %>%
      mutate(log2pkh = log_pkh_imputed/log10(2))
    
    if(m == "interactIBD_otu2"){
      x = "donor_type"
      
      p <- ggplot(data = data, aes(x = as.factor(eval(parse(text = x))), y = log_pkh_imputed)) + 
        geom_violin() +
        geom_boxplot(width = 0.05)  
      
    }
    
    else if (m == "infl_score"){
      xa = m
      
      p <- ggplot(data = data, aes(x = eval(parse(text = xa)),y = log_pkh_imputed)) + 
        geom_bin_2d(bins = 30) +
        geom_smooth(method = "lm")
      
    }
    
    else{
      x = m
      
      p <- ggplot(data = data, aes(x = as.factor(eval(parse(text = x))), y = log_pkh_imputed)) + 
        geom_violin() +
        geom_boxplot(width = 0.05)
      
    }
    
    p +
      labs(title = feat_nm, 
           x = m, y = "Log2 Peak Height",
           subtitle = paste0("Qval: ", 10^-s()$y, "\n", "Avg.Log2FC: ", s()$x)) + 
      theme_bw()
    
  }
  )
  
  # label_name <- eventReactive(s(), {
  #   data.frame(
  #     x = gp()$log2_FC$linmod_df[[2]]$log2_FC,
  #     y = gp()$log2_FC$linmod_df[[2]]$log10_pval_adj,
  #     text = gp()$log2_FC$linmod_df[[2]]$name.y
  #   )[x==s()$log2_FC & y==s()$log10_pval_adj]
  # })
  
  output$feature_violinplot = renderPlot(feat_boxplot())
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)
