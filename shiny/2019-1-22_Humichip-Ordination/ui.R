###############################
#### Humichip Shiny App UI ####
###############################

## Load Libraries -----------------------------------
library(shiny)

# Define UI ----------------------------------------
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Humichip Ordination"),
  
  # Sidebar 
  sidebarLayout(
  sidebarPanel(
      
  # Sample Selection
  fluidRow(column(12,
                  
                  #############
                  ### Visit ###
                  #############
                  checkboxGroupInput('visit_number', 'Visit:', choices = list("Visit 1" = 1,"Visit 4" = 4,"Visit 5" = 5), 
                  selected = c(1, 4, 5), inline = TRUE),
                  
                  ###############
                  ### Matched ###
                  ###############
                  checkboxInput("matched", label = "Matched Samples?", value = TRUE),
                  
                  ########################
                  ### Treatment Groups ###
                  ########################
                  checkboxGroupInput('treatment_groups', 'Treatment Groups:', 
                                     choices = c("LEV", "RIF", "AZI"), 
                                     selected = c("LEV", "RIF", "AZI"), inline = TRUE),
                  
                  ########################
                  ### Detection Method ###
                  ########################
                  br(),
                  h4("Select by Pathogen"),
                  checkboxInput("pathogen_select", label = "Select by Pathogen?", value = FALSE),
                  radioButtons("detection_method", "Detection Method:", 
                               choices = c("taq", "culture", "either", "both"),
                               selected = "both", inline = TRUE),
                  
                  ##########################
                  ### Pathogen Selection ###
                  ##########################
                  checkboxGroupInput("pathogens", "Select Pathogens:",
                                     choices = c("EAEC", "ETEC", "EPEC", "Shigella", "Norovirus"),
                                     selected = c("EAEC", "ETEC", "EPEC", "Shigella", "Norovirus")),
                  checkboxInput("allow_coinfections", label = "Allow Coinfections?", value = FALSE)
  )),
  
  # Perform Ordination Button
  actionButton("action", label = "Analyze"),

  # Side Bar Width
  width = 4),
    
  # Visualize Ordination Plot/data table
  mainPanel(
     fluidRow(column(12,tableOutput("humi_table")))
))))