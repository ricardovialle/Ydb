library(shiny);library(shinyjs);library(markdown);library(DT);library(ggplot2);library(scales)
library(rCharts);library(shinythemes);library(shinyWidgets);library(shinydashboard)
library(shinydashboardPlus);library(shinycssloaders);library(billboarder);library(plotly)
library(doBy);library(plyr);library(dplyr);library(shinydashboardPlus);library(DBI);library(heatmaply)
library(pool);library(datasets);library(reshape2);library(visNetwork);library(d3heatmap);library(ggsci)
library(RMySQL)

home_page = function(){
  fluidPage(
    singleton(tags$head(tags$script(src = "message-handler.js"))),
    column(10, offset = .5,
           h3("Welcome to Ydb! - This page is under construction", style = "font-family: 'Source Sans Pro';"),
           "Choose an option to exhibit the antibody variable domain analysis:",
           h5("Structures from Protein Data Bank (PDB) that contain antibodies:")
    ),
    fluidRow(box(
      title = "Show all stored antibody PDB sequences", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      checkboxGroupInput(inputId = "all_pdb_checkbox", label = "", 
                        choices = list("Show engineered antibodies" = 1, 
                                        "Show only sequences with the right amino acid in conserved positions" = 2,
                                        "Show only one VH/VL pair for each PDB file" = 3,
                                        "Use identity filter" = 4),
                         selected = c(1,2,3)),
      conditionalPanel(
        condition = "input.all_pdb_checkbox.includes('4')",
        sliderInput(inputId = "identity_cutoff", label="% cutoff:", min = 0, max = 100, post  = " %", value = 50, width = '50%')
      ),
      actionButton(inputId = "actionButton.all_pdb_checkbox", label = "Submit"),
      hr(),
      textOutput("all_pdb_checkbox_sel")
    )),
    fluidRow(box(
      title = "Filter by PDB IDs", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      plotOutput("plot1", height = 250)
    )),
    fluidRow(box(
      title = "Filter by antibody-producing organism", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      plotOutput("plot2", height = 250)
    ))
  )
}

#### Header ####
header <- dashboardHeaderPlus(
  title = "Ydb"
)
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
)
#### ####

#### Body ####
body <- dashboardBody(
  tags$head(tags$style(HTML('.main-header .logo { font-size: 16px; }'))),
  tabItems(
    tabItem(tabName = "home",
            home_page()
    )
  )
)
#### ####

  
shinyUI(
  dashboardPagePlus(title = "Ydb",
                    header, 
                    sidebar,
                    body, 
                    collapse_sidebar = T)
)
