###################################
#### Humichip Shiny App Server ####
###################################

## Load Libraries -----------------------------------
library(shiny)
library(tidyverse)
library(vegan)

# Define server -------------------------------------
shinyServer(function(input, output){
  
  # Read in Raw Data with Progress Bar
  withProgress(message = "Reading in Data", min = 0, {
               
               # Humichip
               incProgress(amount = 1/3, detail = "Reading Humichip Data")
               humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip.tsv")))
               
               
               # Treat
               incProgress(amount = 1/3, detail = "Reading Metadata Data")
               treat <- suppressWarnings(suppressMessages(read_csv("TrEAT_Clinical_Metadata_tidy.csv")))
               
               
               # ID Decoder
               incProgress(amount = 1/3, detail = "Reading Decoder")
               ID_decoder <- suppressWarnings(suppressMessages(read_csv("ID_Decoder_Humichip.csv")))
               
               
  })
  
  # Show Humi Data Table
  output$humi_table <- renderTable({head(humichip)})

})
