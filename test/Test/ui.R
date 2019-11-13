library(shiny);library(shinyjs);library(markdown);library(DT);library(ggplot2);library(scales)
library(rCharts);library(shinythemes);library(shinyWidgets);library(shinydashboard)
library(shinydashboardPlus);library(shinycssloaders);library(billboarder);library(plotly)
library(doBy);library(plyr);library(dplyr);library(shinydashboardPlus);library(DBI);library(heatmaply)
library(pool);library(datasets);library(reshape2);library(visNetwork);library(d3heatmap);library(ggsci)
library(RMySQL)

#### Layout for specific pages ####
home_page = function(){
  fluidPage(
    column(10, offset = .5,
           h3("Welcome to Ydb! - This page is under construction", style = "font-family: 'Source Sans Pro';"),
           "Choose an option to exhibit the antibody variable domain analysis:",
           h5("Structures from Protein Data Bank (PDB) that contain antibodies:")
    ), 
    column(10, offset = .5,
           fluidRow(box(
             title = "Show all stored antibody PDB sequences", status = "primary", solidHeader = TRUE,
             collapsible = TRUE, collapsed = T,
             checkboxGroupInput(inputId = "all_pdb_checkbox", label = "", 
                                choices = list("Show engineered antibodies" = 1, 
                                               "Show only sequences with the right amino acid in conserved positions" = 2,
                                               "Show only one VH/VL pair for each PDB file" = 3,
                                               "Use identity filter. Cutoff:" = 4),
                                selected = c(1,2,3)),
             hr(),
             textOutput("all_pdb_checkbox_sel")
           ))
           
    )
  )
}

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  #### Dashboard page layout parts ----
  # Include: Header + Sidebars + Body
  
  #### Header ####
  header <- dashboardHeaderPlus(
    title = "Ydb"
  ),
  #### ####
  
  #### Sidebar ####
  sidebar <- dashboardSidebar(
    tags$head(tags$style(type = "text/css", ".shiny-input-container {padding-top: 0px !important;padding-bottom: 0px !important;}")),
    
    sidebarMenu(id = "menu",
                menuItem("Home", tabName = "home", icon = icon("home")),
                menuItem("Database explorer", tabName = "database", icon = icon("database"),
                         menuSubItem("by Complex", tabName = "complex"),
                         menuSubItem("by Chain", tabName = "chain"),
                         menuSubItem("by Interaction", tabName = "interaction")
                ),
                menuItem("Visualization", tabName = "charts_menu", icon = icon("bar-chart-o"),
                         menuSubItem("Database Statistics", tabName = "dbstats_charts"),
                         menuSubItem("Interactions", tabName = "interaction_charts"),
                         menuSubItem("Custom Graphs", tabName = "custom_charts")
                )
    )
  ),
  #### ####
  
  #### Body ####
  body <- dashboardBody(
    tags$head(tags$style(HTML('.main-header .logo { font-size: 16px; }'))),
    tabItems(
      tabItem(tabName = "home",
              home_page()
      )
    )
  ),
  #### ####
  
  #### Page call (main function) ----
  dashboardPagePlus(title = "Ydb",
                    header, 
                    sidebar, 
                    body, 
                    collapse_sidebar = T)
  #### ####
  
))
