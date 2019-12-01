library(shiny);library(shinyjs);library(markdown);library(DT);library(ggplot2);library(scales)
library(rCharts);library(shinythemes);library(shinyWidgets);library(shinydashboard)
library(shinydashboardPlus);library(shinycssloaders);library(billboarder);library(plotly)
library(doBy);library(plyr);library(dplyr);library(shinydashboardPlus);library(DBI);library(heatmaply)
library(pool);library(datasets);library(reshape2);library(visNetwork);library(d3heatmap);library(ggsci)
library(RMySQL)

#################################################################################
#### Definitions ################################################################
#################################################################################

#### Atom interaction ####
atom_interaction_type = data.frame(
  code = c("F","E","A","P","V","C","H","I","J","K","L","M","N","O","Q","R","S","T"),
  type = c("Hydrophobic","Electrostatic","Aromatic","Pi-cation",
           "Van der Waals","Van der Waals clash","Hydrogen bond-MM",
           "Hydrogen bond-SM","Hydrogen bond-MS","Hydrogen bond-SS",
           "Hydrogen bond-Water mediated-MM","Hydrogen bond-Water mediated-SM",
           "Hydrogen bond-Water mediated-MS","Hydrogen bond-Water mediated-SS",
           "Hydrogen bond-2 Waters mediated-MM","Hydrogen bond-2 Waters mediated-SM",
           "Hydrogen bond-2 Waters mediated-MS","Hydrogen bond-2 Waters mediated-SS")
)
#### ####

#### AA classes #####
aa_classes = data.frame(
  class = c("Acid","Acid",
            "Basic","Basic","Basic",
            "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
            "Neutral","Neutral",
            "Polar","Polar","Polar","Polar","Polar",
            "Gap"
  ), 
  class_color = c("red","red",
                  "blue","blue","blue",
                  "black","black","black","black","black","black","black","black",
                  "purple","purple",
                  "green","green","green","green","green",
                  "gray"
  ),
  aa = c("D","E",
         "R","H","K",
         "F","P","V","L","I","A","W","M",
         "N","Q",
         "C","S","G","Y","T",
         "-")
)
#### ####

#### Complex table - Columns dictionary ####
complex_tbl_dic = c(
  "Id","Complex Id","Complex",
  "Pdb","PDB","Complex",
  "Model","Model","Complex",
  "DeltaG","Delta G","Complex",
  "DeltaSAS","Delta SAS","Complex",
  "Affinity","Affinity","Complex",
  "AffinityMethod","Affinity Method","Complex",
  "AffinityTemperature","Affinity Temperature","Complex",
  "AffinityPmid","Affinity PMID","Complex",
  "PDBDescription","PDB Description","Complex",
  "PDBCompound","PDB Compound","Complex",
  "PDBAuthors","PDB Authors","Complex",
  "PDBReference","PDB Reference","Complex",
  "PDBDepositionDate","PDB Deposition Date","Complex",
  "PDBExperiment","PDB Experiment","Complex",
  "PDBResolution","PDB Resolution","Complex",
  "PDBRFree","R Free","Complex",
  "PDBRFactor","R Factor","Complex",
  "HeavyChainId","Chain Id (Heavy)","Chain",
  "HeavyChainName","Chain Name (Heavy)","Chain",
  "HeavyChainType","Chain Type (Heavy)","Chain",
  "HeavyChainSubclass","Chain Subclass (Heavy)","Chain",
  "HeavySpecies","Chain Species (Heavy)","Chain",
  "HeavyLightType","Light Type (Heavy)","Chain",
  "HeavyResSequence","Sequence (Heavy)","Chain",
  "HeavyAtomSequence","Atom Sequence (Heavy)","Chain",
  "HeavyGravy","Gravy (Heavy)","Chain",
  "HeavyIsoelectricPoint","Isoelectric Point (Heavy)","Chain",
  "HeavyNumberOfAas","Number of Amino Acids (Heavy)","Chain",
  "HeavyTotalNumberofAas","Total Number of AAs (Heavy)","Chain",
  "HeavyAliphaticIndex","Aliphatic Index (Heavy)","Chain",
  "HeavyMolecularWeight","Molecular Weight (Heavy)","Chain",
  "HeavyNegativeCharged","Negative Charged (Heavy)","Chain",
  "HeavyPositiveCharged","Positive Charged (Heavy)","Chain",
  "HeavyPolarUncharged","Polar Uncharged (Heavy)","Chain",
  "HeavySmall","Small (Heavy)","Chain",
  "HeavyHydrophobic","Hydrophobic (Heavy)","Chain",
  "HeavyAlcohol","Alcohol (Heavy)","Chain",
  "HeavyAromatic","Aromatic (Heavy)","Chain",
  "HeavyAlphaHelixStructure","Alpha-Helix (Heavy)","SS",
  "HeavyIsolatedBetaBridgeStructure","Isolated Beta Bridge (Heavy)","SS",
  "HeavyStrandStructure","Strand (Heavy)","SS",
  "HeavyHelix3_10Structure","Helix 3_10 (Heavy)","SS",
  "HeavyPiHelixStructure","Pi Helix (Heavy)","SS",
  "HeavyTurnStructure","Turn (Heavy)","SS",
  "HeavyBendStructure","Bend (Heavy)","SS",
  "HeavyNoneStructure","None secondary structure (Heavy)","SS",
  "HeavyIMGTComment","IMGT Comment (Heavy)","Chain",
  "HeavyPdbComment","PDB Comment (Heavy)","Chain",
  "LightChainId","Chain Id (Light)","Chain",
  "LightChainName","Chain Name (Light)","Chain",
  "LightChainType","Chain Type (Light)","Chain",
  "LightChainSubclass","Chain Subclass (Light)","Chain",
  "LightSpecies","Chain Species (Light)","Chain",
  "LightLightType","Light Type (Light)","Chain",
  "LightResSequence","Sequence (Light)","Chain",
  "LightAtomSequence","Atom Sequence (Light)","Chain",
  "LightGravy","Gravy (Light)","Chain",
  "LightIsoelectricPoint","Isoelectric Point (Light)","Chain",
  "LightNumberOfAas","Number of Amino Acids (Light)","Chain",
  "LightTotalNumberofAas","Total Number of AAs (Light)","Chain",
  "LightAliphaticIndex","Aliphatic Index (Light)","Chain",
  "LightMolecularWeight","Molecular Weight (Light)","Chain",
  "LightNegativeCharged","Negative Charged (Light)","Chain",
  "LightPositiveCharged","Positive Charged (Light)","Chain",
  "LightPolarUncharged","Polar Uncharged (Light)","Chain",
  "LightSmall","Small (Light)","Chain",
  "LightHydrophobic","Hydrophobic (Light)","Chain",
  "LightAlcohol","Alcohol (Light)","Chain",
  "LightAromatic","Aromatic (Light)","Chain",
  "LightAlphaHelixStructure","Alpha-Helix (Light)","SS",
  "LightIsolatedBetaBridgeStructure","Isolated Beta Bridge (Light)","SS",
  "LightStrandStructure","Strand (Light)","SS",
  "LightHelix3_10Structure","Helix 3_10 (Light)","SS",
  "LightPiHelixStructure","Pi Helix (Light)","SS",
  "LightTurnStructure","Turn (Light)","SS",
  "LightBendStructure","Bend (Light)","SS",
  "LightNoneStructure","None secondary structure (Light)","SS",
  "LightIMGTComment","IMGT Comment (Light)","Chain",
  "LightPdbComment","PDB Comment (Light)","Chain",
  "AntigenName","Antigen Name","Chain",
  "AntigenChainId","Chain Id (Antigen)","Chain",
  "AntigenChainName","Chain Name (Antigen)","Chain",
  "AntigenChainType","Chain Type (Antigen)","Chain",
  "AntigenChainSubclass","Chain Subclass (Antigen)","Chain",
  "AntigenSpecies","Chain Species (Antigen)","Chain",
  "AntigenLightType","Light Type (Antigen)","Chain",
  "AntigenResSequence","Sequence (Antigen)","Chain",
  "AntigenAtomSequence","Atom Sequence (Antigen)","Chain",
  "AntigenGravy","Gravy (Antigen)","Chain",
  "AntigenIsoelectricPoint","Isoelectric Point (Antigen)","Chain",
  "AntigenNumberOfAas","Number of Amino Acids (Antigen)","Chain",
  "AntigenTotalNumberofAas","Total Number of AAs (Antigen)","Chain",
  "AntigenAliphaticIndex","Aliphatic Index (Antigen)","Chain",
  "AntigenMolecularWeight","Molecular Weight (Antigen)","Chain",
  "AntigenNegativeCharged","Negative Charged (Antigen)","Chain",
  "AntigenPositiveCharged","Positive Charged (Antigen)","Chain",
  "AntigenPolarUncharged","Polar Uncharged (Antigen)","Chain",
  "AntigenSmall","Small (Antigen)","Chain",
  "AntigenHydrophobic","Hydrophobic (Antigen)","Chain",
  "AntigenAlcohol","Alcohol (Antigen)","Chain",
  "AntigenAromatic","Aromatic (Antigen)","Chain",
  "AntigenAlphaHelixStructure","Alpha-Helix (Antigen)","SS",
  "AntigenIsolatedBetaBridgeStructure","Isolated Beta Bridge (Antigen)","SS",
  "AntigenStrandStructure","Strand (Antigen)","SS",
  "AntigenHelix3_10Structure","Helix 3_10 (Antigen)","SS",
  "AntigenPiHelixStructure","Pi Helix (Antigen)","SS",
  "AntigenTurnStructure","Turn (Antigen)","SS",
  "AntigenBendStructure","Bend (Antigen)","SS",
  "AntigenNoneStructure","None secondary structure (Antigen)","SS",
  "AntigenIMGTComment","IMGT Comment (Antigen)","Chain",
  "AntigenPdbComment","PDB Comment (Antigen)","Chain",
  "Interaction","Interaction","Complex"
)
complex_tbl_dic = data.frame(colname = complex_tbl_dic[seq(1,length(complex_tbl_dic),by = 3)],
                             desc = complex_tbl_dic[seq(2,length(complex_tbl_dic),by = 3)],
                             class = complex_tbl_dic[seq(3,length(complex_tbl_dic),by = 3)])
#### ####

#### Chain table - Columns dictionary ####
chain_tbl_dic = c(
  "Id","Chain Id","Chain",
  "Pdb","PDB","Chain",
  "ChainName","Chain Name","Chain",
  "ChainType","Chain Type","Chain",
  "SpeciesId","Species Id","Chain",
  "Species","Species","Chain",
  "Gravy","Gravy","Chain",
  "IsoelectricPoint","Isoelectric Point","Chain",
  "AliphaticIndex","Aliphatic Index","Chain",
  "MolecularWeight","Molecular Weight","Chain",
  "TotalNumberOfAas","Number of AAs (Total)","Chain",
  "NumberOfAas","Number of AAs","Chain",
  "NegativeCharged","Negative Charged","Chain",
  "PositiveCharged","Positive Charged","Chain",
  "PolarUncharged","Polar Uncharged","Chain",
  "Small","Small","Chain",
  "Hydrophobic","Hydrophobic","Chain",
  "Alcohol","Alcohol","Chain",
  "Aromatic","Aromatic","Chain",
  "AlphaHelixStructure","Alpha-Helix","SS",        
  "IsolatedBetaBridgeStructure","Isolated Beta Bridge","SS",
  "StrandStructure","Strand","SS",
  "Helix3_10Structure","Helix 3-10","SS",
  "PiHelixStructure","Pi Helix","SS",
  "TurnStructure","Turn","SS",
  "BendStructure","Bend","SS",
  "NoneStructure","None secondary structure","SS"        
)
chain_tbl_dic = data.frame(colname = chain_tbl_dic[seq(1,length(chain_tbl_dic),by = 3)],
                           desc = chain_tbl_dic[seq(2,length(chain_tbl_dic),by = 3)],
                           class = chain_tbl_dic[seq(3,length(chain_tbl_dic),by = 3)])
#### ####

#### Interaction table - Columns dictionary ####
interaction_tbl_dic = c(
  "AtomInteraction_AntibodyAtom","Antibody Atom","Antibody",
  "AtomInteraction_AntibodyResidueId","Antibody Residue","Antibody",
  "AntibodyResidue_Id","Antibody Res Id","Antibody",
  "AntibodyResidue_ChainId","Antibody Chain Id","Antibody",
  "AntibodyResidue_PdbIndex","Residue PDB Index (Antibody)","Antibody",
  "AntibodyResidue_PdbIndexCode","Residue PDB Index code (Antibody)","Antibody",
  "AntibodyResidue_AminoAcid","Amino acid (Antibody)","Antibody",
  "AntibodyResidue_Flexibility","Residue Flexibility (Antibody)","Antibody",
  "AntibodyResidue_FlexibilityConfidence","Flexibility Confidence (Antibody)","Antibody",
  "AntibodyResidue_StructuralProperty","Residue Structural Property (Antibody)","Antibody",
  "AntibodyResidue_SASComplexed","Residue SAS Complexed (Antibody)","Antibody",
  "AntibodyResidue_SASUncomplexed","Residue SAS Uncomplexed (Antibody)","Antibody",
  "AntibodyResidue_IMGTIndex","Residue IMGT Index (Antibody)","Antibody",
  "AtomInteraction_AntigenAtom","Antigen Atom","Antigen",
  "AtomInteraction_AntigenResidueId","Antigen Residue","Antigen",
  "AntigenResidue_Id","Antigen Res Id","Antigen",
  "AntigenResidue_ChainId","Antigen Chain Id","Antigen",
  "AntigenResidue_PdbIndex","Residue PDB Index (Antigen)","Antigen",
  "AntigenResidue_PdbIndexCode","Residue PDB Index code (Antigen)","Antigen",
  "AntigenResidue_AminoAcid","Amino acid (Antigen)","Antigen",
  "AntigenResidue_Flexibility","Residue Flexibility (Antigen)","Antigen",
  "AntigenResidue_FlexibilityConfidence","Flexibility Confidence (Antigen)","Antigen",
  "AntigenResidue_StructuralProperty","Residue Structural Property (Antigen)","Antigen",
  "AntigenResidue_SASComplexed","Residue SAS Complexed (Antigen)","Antigen",
  "AntigenResidue_SASUncomplexed","Residue SAS Uncomplexed (Antigen)","Antigen",
  "AntigenResidue_IMGTIndex","Residue IMGT Index (Antigen)","Antigen",
  "AtomInteraction_Type","Interaction Type","Interaction",
  "AtomInteraction_Distance","Interaction Distance","Interaction"
)
interaction_tbl_dic = data.frame(colname = interaction_tbl_dic[seq(1,length(interaction_tbl_dic),by = 3)],
                                 desc = interaction_tbl_dic[seq(2,length(interaction_tbl_dic),by = 3)],
                                 class = interaction_tbl_dic[seq(3,length(interaction_tbl_dic),by = 3)])
#### ####

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
