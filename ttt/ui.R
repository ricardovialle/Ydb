library(shiny);library(shinyjs);library(markdown);library(DT);library(ggplot2);library(scales)
library(rCharts);library(shinythemes);library(shinyWidgets);library(shinydashboard)
library(shinydashboardPlus);library(shinycssloaders);library(billboarder);library(plotly)
library(doBy);library(plyr);library(dplyr);library(shinydashboardPlus);library(DBI);library(heatmaply)
library(pool);library(datasets);library(reshape2);library(visNetwork);library(d3heatmap);library(ggsci)
library(RMySQL)

#### Home page ####
home_page <- function(){
  fluidPage(
    singleton(tags$head(tags$script(src = "message-handler.js"))),
    column(10, offset = .5,
           h3("Welcome to Ydb! - This page is under construction", style = "font-family: 'Source Sans Pro';"),
           "Choose an option to exhibit the antibody variable domain analysis:",
           h5("Structures from Protein Data Bank (PDB) that contain antibodies:")
    ),
    
    #### Show all box ####
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
      hr()
    )),
    #### ####
    
    #### Filter by PDB ####
    fluidRow(box(
      title = "Filter by PDB IDs", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      radioButtons(inputId = "filter_by_pdb_radiobutton", label = "", inline = T, 
                         choices = list("Specify PDB IDs" = 1, 
                                        "Specify PDB IDs and chain name" = 2)),
      conditionalPanel(
        condition = "input.filter_by_pdb_radiobutton=='1'",
        textAreaInput(inputId = "filter_by_pdb_textarea", label="Enter PDB IDs separated by commas, semicolons or in new lines:", value = "", width = NULL,
                      height = NULL, cols = NULL, rows = NULL, placeholder = NULL,
                      resize = NULL)
      ),
      conditionalPanel(
        condition = "input.filter_by_pdb_radiobutton=='2'",
        textAreaInput(inputId = "filter_by_pdb_textarea", label="Enter a list of chains separated by commas, semicolons or in new lines. Each chain must be specified by PDB ID followed by colon and chain name:", value = "", width = NULL,
                      height = NULL, cols = NULL, rows = NULL, placeholder = NULL,
                      resize = NULL)
      ),
      actionButton(inputId = "actionButton.filter_by_pdb_radiobutton", label = "Submit"),
      hr()
    )),
    #### ####
    
    #### Filter by antibody organism ####
    fluidRow(box(
      title = "Filter by antibody-producing organism", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      
      selectInput(
        "filter_by_antibody_organism_selectinput", "Choose one or more species from the list:", 
        c("Human","Cow","Dog"),
        multiple = TRUE
      ),
      actionButton(inputId = "actionButton.filter_by_antibody_organism_selectinput", label = "Submit"),
      hr()
    )),
    #### ####
    
    #### Filter by antigen organism ####
    fluidRow(box(
      title = "Filter by antigen-producing organism", status = "primary", solidHeader = TRUE,
      collapsible = TRUE, collapsed = T,
      
      selectInput(
        "filter_by_antigen_organism_selectinput", "Choose one or more species from the list:", 
        c("Human","Cow","Dog"),
        multiple = TRUE
      ),
      actionButton(inputId = "actionButton.filter_by_antigen_organism_selectinput", label = "Submit"),
      hr()
    ))
    #### ####
    
  )
}
#### ####

#### Main Tables page ####
tables_page <- function(){
  fluidPage(
    
    #### Opts ####
    theme = shinytheme("yeti"),
    useShinyjs(),
    #### ####
    
    singleton(tags$head(tags$script(src = "message-handler.js"))),
    column(10, offset = .5,
           h3("Welcome to Ydb! - This page is under construction", style = "font-family: 'Source Sans Pro';"),
           "Select how you would like to view:",
           h5("##### SEARCH PARAMETERS #####"),
           textOutput("all_pdb_checkbox_sel"),
           br()
    ),
    #### Row on top - (options) ####
    fixedRow(
      column(12            
             ,div(style="display:inline-block",downloadButton("downloadData", "Download",class = "btn btn-primary"))
             ,div(style="display:inline-block",dropdownButton(inputId = "columnsToShow", label = "Show columns",
                                                              circle = F, status = "primary",
                                                              selectizeInput("complex_tbl_columnSelection_by_complex", label = h5("PDB"), multiple = T, width = '100%',
                                                                             choices = as.character(complex_tbl_dic[which(complex_tbl_dic$class=="Complex"),c("desc")]),
                                                                             selected = "PDB"),
                                                              selectizeInput("complex_tbl_columnSelection_by_chain", label = h5("Chain"), multiple = T, 
                                                                             choices = as.character(complex_tbl_dic[which(complex_tbl_dic$class=="Chain"),c("desc")]),
                                                                             selected = c("Chain Name (Heavy)","Chain Name (Light)","Chain Name (Antigen)")),
                                                              selectizeInput("complex_tbl_columnSelection_by_ss", label = h5("Secondary Structure"), multiple = T,
                                                                             choices = as.character(complex_tbl_dic[which(complex_tbl_dic$class=="SS"),c("desc")]))
             ))
             ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
             ,div(style="display:inline-block",actionButton("showComplexesForSelection","See complexes",class = "btn btn-success"))
             ,div(style="display:inline-block",actionButton("showChainsForSelection","See chains",class = "btn btn-success"))
             ,div(style="display:inline-block",actionButton("showInteractionsForSelection","See interactions",class = "btn btn-success"))
             ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
             ,div(style="display:inline-block",dropdownButton(inputId = "complex_table_seegraphs", label = "Plot interactions",
                                                              circle = F, status = "warning", icon = icon("bar-chart-o"), 
                                                              actionButton("select_antibodyresidue_by_position_barplot","Frequency of Antibody amino acids by IMGT position",class = "btn btn-warning",width = '100%'),
                                                              actionButton("select_atominteraction_residues_heatmap","Number of interactions by antibody-antigen residues",class = "btn btn-warning", width = '100%'),
                                                              actionButton("select_atominteraction_residues_network","Network",class = "btn btn-warning", width = '100%')))
             ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
             ,div(style="display:inline-block;text-align:right",dropdownButton(
               tags$h4("Selection info"),
               div(style="font-size: small;font-weight: normal",htmlOutput("selection_info_complexTableView")),hr(),
               div(style="font-size: small;font-weight: normal",materialSwitch(inputId = "complex_selection_switch",label = "Show only selected entries", status = "primary", right = F)),
               circle = TRUE, status = "info", icon = icon("info")
             ))
      )
    )
    #### ####
    ,br(),
    mainPanel(
      tabsetPanel(
        tabPanel("by Complex", DT::dataTableOutput("complex_table")),
        tabPanel("by Chain", DT::dataTableOutput("mytable2")),
        tabPanel("by Interaction", DT::dataTableOutput("mytable3"))
      )
    )
  )
}
#### ####

#### Stats plots page ####
dbstats_charts = function(){
  fluidPage(
    #### opts ####
    theme = shinytheme("yeti"),
    useShinyjs(),
    
    box(
      title = "General Statistics",
      status = "warning",
      width = 12,
      footer = withSpinner(tableOutput("general_stats_tbl"))
    ),
    
    box(
      title = "PDB Statistics",
      status = "info",
      width = 12,
      collapsible = T, 
      collapsed = F,
      billboarderOutput(outputId = "pdb_byyear_barplot"),
      billboarderOutput(outputId = "pdbexperiment_piechart"),
      billboarderOutput(outputId = "pdb_resolution_histogram")
    ),
    
    box(
      title = "Chain Statistics",
      status = "info",
      width = 12,
      collapsible = T, 
      collapsed = F,
      billboarderOutput(outputId = "heavy_species_piechart"),
      billboarderOutput(outputId = "light_species_piechart"),
      billboarderOutput(outputId = "antigen_species_piechart")
    )
    
  )
}
#### ####

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
    ),
    tabItem(tabName = "complex",
            tables_page()
    ),
    tabItem(tabName = "dbstats_charts",
            dbstats_charts()
    )
  )
)
#### ####

#### Page call (main function) ----
shinyUI(
  dashboardPagePlus(title = "Ydb",
                    header, 
                    sidebar,
                    body, 
                    collapse_sidebar = T)
)
