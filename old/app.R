#### Libraries / Database connection / Global variables ###########################
source("global.R") # Run once as app start 
###################################################################################

# Suppress warning messages on console
options(shiny.sanitize.errors = FALSE, warn =-1)

# Global variables
complex_table = complex_table_show_all()
PdbList = as.factor(sort(unique(complex_table$Pdb)))
ExperimentList = as.factor(sort(unique(complex_table$PDBExperiment)))
SpeciesList = as.factor(sort(unique(c(complex_table$HeavySpecies,complex_table$LightSpecies,complex_table$AntigenSpecies))))
chain_table = chain_table_show_all()
interaction_table = interaction_table_show_all()

###################################################################################
# Server logic ----
server <- function(input, output, session) {

  # Automatically stop a Shiny app when closing the browser tab
  #session$onSessionEnded(stopApp)
  
  # Set this to "force" instead of TRUE for testing locally (without Shiny Server)
  #session$allowReconnect("force")
  #session$allowReconnect(TRUE)
  
  # Default definitions
  columns_to_show_complex_table = c("Pdb","HeavyChainName", "LightChainName", "AntigenChainName")
  columns_to_show_chain_table = c("Pdb","ChainName","ChainType","Species")
  columns_to_show_interaction_table = c("AtomInteraction_AntibodyAtom","AntibodyResidue_AminoAcid",
                                        "AtomInteraction_AntigenAtom","AntigenResidue_AminoAcid",
                                        "AtomInteraction_Type","AtomInteraction_Distance")
  
  #################################################################################
  #### HOME: Selection Menus ######################################################
  #################################################################################

  output$all_pdb_checkbox_sel <- renderPrint({ input$all_pdb_checkbox })

    
  #################################################################################
  #### Selection Control ##########################################################
  ################################################################################# 
  # All tables (Complex, Chain and Interaction) can be connected in the background
  # This permits selection and filtering across each other. 
  # The following codes are to control it programmatically.  
  
  #### Reactive values to store selected rows ####
  # This is a reactiveValue: When you read a value from it, 
  # the calling reactive expression takes a reactive dependency on that value, 
  # and when you write to it,
  # it notifies any reactive functions that depend on that value.
  
  # For complex_table
  selection_all <- reactiveValues( 
    all = data.frame( # This contains a data.frame with all rownames in table ('ID') and which one is selected ('sel')
      ID = rownames(complex_table), # Using rownames as reference *depends on the MySQL select' output order
      sel = rep(F,nrow(complex_table))
    )
  )

  # For chain_table (same as above)
  selection_all_chain <- reactiveValues(
    all = data.frame(
      ID = chain_table$Id, # Using the ChainID as reference
      sel = rep(F,nrow(chain_table))
    )
  )
  
  # interaction_table does not have this control
  # it relies on the selection in the complex/chain tables
  
  #### ####
  
  #### Currently selected dataset ####
  curSelectedData_complex_table <- reactive({
    # See what is selected
    rows_selected = which(selection_all$all$sel==T)
    if (length(rows_selected)==0){ # If nothing has been selected show all entries
      table_selected = complex_table
    }else{ # If something has been selected, show only selection 
      table_selected = complex_table[rows_selected,]
    }
    #switch(input$dataset1, All = complex_table, Selected = table_selected)
  })
  #### ####
  
  #### Currently selected dataset ####
  curSelectedData_chain_table <- reactive({
    # See what is selected
    rows_selected = which(selection_all_chain$all$sel==T)
    if (length(rows_selected)==0){ # If nothing has been selected show all entries
      table_selected = chain_table
    }else{ # If something has been selected, show only selection 
      table_selected = chain_table[rows_selected,]
    }
    #switch(input$dataset1, All = chain_table, Selected = table_selected)
  })
  #### ####
  
  #### Observer for selected rows ####
  # For complex_table
  observeEvent(input$tbl_cell_clicked$row,{
    # This controls case a diferent data is being shown (e.g. after filtering)
    data = curFilteredData_complex_table()
    a_sel <- isolate(selection_all$all$sel[selection_all$all$ID %in% rownames(data)]) # Get what is already selected given the current data
    a <- data.frame(IDs = rownames(data), sel = a_sel) # Temporary data.frame with selections for the current data
    # Control of click (selection) on the current table
    if (is.null(input$tbl_cell_clicked$row)) {
      a[, 'sel'] <- F
      return()
    } 
    if (isTRUE(a[input$tbl_cell_clicked$row, 'sel'])){ # If complex was already selected. Deselect it.
      a[input$tbl_cell_clicked$row, 'sel'] <- F
      # Update values in the complete version of table
      selection_all$all$sel[selection_all$all$ID %in% a$IDs] <- a$sel
      # Update values for chains and interactions as well...
      complex_rows_deselected = a[input$tbl_cell_clicked$row, 'IDs'] # Deselected complex
      chain_ids = complex_table[which(rownames(complex_table)==complex_rows_deselected),c("HeavyChainId","LightChainId","AntigenChainId")] # Look chain ids in the complex table
      selection_all_chain$all$sel[selection_all_chain$all$ID %in% chain_ids] <- F # Deselected chains related to complexes
      return()
    }
    if (!isTRUE(a[input$tbl_cell_clicked$row, 'sel'])){
      a[input$tbl_cell_clicked$row, 'sel'] <- T
      # Update values in the complete version of table
      selection_all$all$sel[selection_all$all$ID %in% a$IDs] <- a$sel
      # Update values for chains and interactions as well...
      complex_rows_selected = isolate(which(selection_all$all$sel==T)) # Get a list of selected complexes
      chain_ids = complex_table[complex_rows_selected,c("HeavyChainId","LightChainId","AntigenChainId")] # Look chain ids in the complex table
      chain_ids = chain_ids[!is.na(chain_ids)] # Remove NA entries
      selection_all_chain$all$sel[selection_all_chain$all$ID %in% chain_ids] <- T # Set selected for chains related to complexes
      return()
    }
  }, ignoreNULL = TRUE)
  # For chain_table
  observeEvent(input$chain_tbl_cell_clicked$row,{
    # This controls case a diferent data is being shown (e.g. after filtering)
    data = curFilteredData_chain_table()
    a_sel <- isolate(selection_all_chain$all$sel[selection_all_chain$all$ID %in% data$Id])
    a <- data.frame(IDs = data$Id, sel = a_sel)
    if (is.null(input$chain_tbl_cell_clicked$row)) {
      a[, 'sel'] <- F
    } else if (isTRUE(a[input$chain_tbl_cell_clicked$row, 'sel'])){
      a[input$chain_tbl_cell_clicked$row, 'sel'] <- F
    } else  if (!isTRUE(a[input$chain_tbl_cell_clicked$row, 'sel'])){
      a[input$chain_tbl_cell_clicked$row, 'sel'] <- T
    }
    selection_all_chain$all$sel[selection_all_chain$all$ID %in% a$IDs] <- a$sel
  }, ignoreNULL = TRUE)
  #### ####
  
  #### To be removed : Notification - Selected ####
  #output$notif <- renderMenu({
  #  ids = which(selection_all$all$sel==T)
  #  nids = length(ids)
  #  chain_ids = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
  #  nchains = length(chain_ids)
  #  # If some chain was selected, count interactions
  #  if (nchains==0){
  #    ninteractions = 0
  #  }else{
  #    data <- interaction_table_by_chainid(chain_ids)
  #    ninteractions = nrow(data)
  #  }
  #  # Notification
  #  if (nids==1){
  #    dropdownMenu(type = "notifications", 
  #                 icon = icon("shopping-cart", lib = "glyphicon"),
  #                 notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                  paste(nids,'complex selected')
  #                 ),
  #                 notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                  paste(nchains,'chains selected')
  #                 ),
  #                 notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                  paste(ninteractions,'interactions selected')
  #                 )
  #    )
  #  }else{
  #    if (nids>1){
  #      dropdownMenu(type = "notifications", 
  #                   icon = icon("shopping-cart", lib = "glyphicon"),
  #                   notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                    paste(nids,'complexes selected')
  #                   ),
  #                   notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                    paste(nchains,'chains selected')
  #                   ),
  #                   notificationItem(icon = icon("ok", lib = "glyphicon"), status = "danger",
  #                                    paste(ninteractions,'interactions selected')
  #                   )
  #      )
  #    }else{
  #      dropdownMenu(type = "notifications", 
  #                   icon = icon("shopping-cart", lib = "glyphicon"))
  #    }
  #  }
  #})
  #### ####
  
  #### observer of select all button ####
  # complex_table
  observeEvent(input$select_all_button, {
    data <- curFilteredData_complex_table()
    a_sel <- selection_all$all$sel[selection_all$all$ID %in% rownames(data)]
    a <- data.frame(IDs = rownames(data), sel = a_sel)
    a[input$tbl_rows_all, 'sel'] <- T
    selection_all$all$sel[selection_all$all$ID %in% a$IDs] <- a$sel
    # Update values for chains and interactions as well...
    complex_rows_selected = isolate(which(selection_all$all$sel==T)) # Get a list of selected complexes
    chain_ids = complex_table[complex_rows_selected,c("HeavyChainId","LightChainId","AntigenChainId")] # Look chain ids in the complex table
    chain_ids = chain_ids[!is.na(chain_ids)] # Remove NA entries
    selection_all_chain$all$sel[selection_all_chain$all$ID %in% chain_ids] <- T # Set selected for chains related to complexes
  })
  # chain_table
  observeEvent(input$select_all_button_chain, {
    data <- curFilteredData_chain_table()
    a_sel <- selection_all_chain$all$sel[selection_all_chain$all$ID %in% data$Id]
    a <- data.frame(IDs = data$Id, sel = a_sel)
    a[input$chain_tbl_rows_all, 'sel'] <- T
    selection_all_chain$all$sel[selection_all_chain$all$ID %in% a$IDs] <- a$sel
  })
  # interaction_table
  observeEvent(input$select_all_button_interaction, {
    data <- curFilteredData_interaction_table()
    a_sel <- selection_all_interaction$all$sel[selection_all_interaction$all$ID %in% rownames(data)]
    a <- data.frame(IDs = rownames(data), sel = a_sel)
    a[input$interaction_tbl_rows_all, 'sel'] <- T
    selection_all_interaction$all$sel[selection_all_interaction$all$ID %in% a$IDs] <- a$sel
  })
  #### ####
  
  #### observer of deselect all button ####
  # complex_table
  observeEvent(input$deselect_all_button, {
    data <- curFilteredData_complex_table()
    a_sel <- selection_all$all$sel[selection_all$all$ID %in% rownames(data)]
    a <- data.frame(IDs = rownames(data), sel = a_sel)
    a[input$tbl_rows_all, 'sel'] <- F
    selection_all$all$sel[selection_all$all$ID %in% a$IDs] <- a$sel
    # Update values for chains and interactions as well...
    complex_rows_deselected = a[input$tbl_cell_clicked$row, 'IDs'] # Deselected complex
    chain_ids = complex_table[which(rownames(complex_table)==complex_rows_deselected),c("HeavyChainId","LightChainId","AntigenChainId")] # Look chain ids in the complex table
    selection_all_chain$all$sel[selection_all_chain$all$ID %in% chain_ids] <- F # Deselected chains related to complexes
  })
  # chain_table
  observeEvent(input$deselect_all_button_chain, {
    data <- curFilteredData_chain_table() 
    a_sel <- selection_all_chain$all$sel[selection_all_chain$all$ID %in% data$Id]
    a <- data.frame(IDs = data$Id, sel = a_sel)
    a[input$chain_tbl_rows_all, 'sel'] <- F
    selection_all_chain$all$sel[selection_all_chain$all$ID %in% a$IDs] <- a$sel
  })
  # interaction_table
  #observeEvent(input$deselect_all_button_interaction, {
  #  data <- curFilteredData_interaction_table()
  #  a_sel <- selection_all_interaction$all$sel[selection_all_interaction$all$ID %in% rownames(data)]
  #  a <- data.frame(IDs = rownames(data), sel = a_sel)
  #  a[input$interaction_tbl_rows_all, 'sel'] <- F
  #  selection_all_interaction$all$sel[selection_all_interaction$all$ID %in% a$IDs] <- a$sel
  #})
  #### ####
  
  #################################################################################
  #### Filtering Control ##########################################################
  ################################################################################# 
  # Each table webpage has it own filtering options
  
  #### Current dataset (controls data after filtering) ####
  # For complex_table
  curFilteredData_complex_table <- reactive({
    data = complex_table # Table is already loaded (from global.R)
    # Filter by PDB
    if (is.null(input$filter_pdb)){
      data <- data
    }  else if (any(input$filter_pdb != "")){
      data <- data[which(data$Pdb %in% input$filter_pdb),]
    }
    # Filter by PDB Experiment
    if (is.null(input$filter_experiment)){
      data <- data
    }  else if (any(input$filter_experiment != "")){
      data <- data[which(data$PDBExperiment%in%input$filter_experiment),]
    }
    # Filter by PDB Resolution
    if (is.null(input$filter_resolution)){
      data <- data
    }  else if (any(input$filter_resolution != "")){
      min_range = min(input$filter_resolution)
      max_range = max(input$filter_resolution)
      filt_idx = which(data$PDBResolution >= min_range & data$PDBResolution <= max_range)
      data = data[filt_idx,]
    }
    # Filter by Species
    if (is.null(input$filter_species)){
      data <- data
    }  else if (any(input$filter_species != "")){
      data = data[which(data$HeavySpecies%in%input$filter_species | data$LightSpecies%in%input$filter_species | data$AntigenSpecies%in%input$filter_species),]
    }
    return(data)
  })
  # For chain_table
  curFilteredData_chain_table <- reactive({
    data = chain_table # Table is already loaded (from global.R)
    # Filter by PDB
    if (is.null(input$chain_filter_pdb)){
      data <- data
    }  else if (any(input$chain_filter_pdb != "")){
      data <- data[which(data$Pdb%in%input$chain_filter_pdb),]
    }
    # Filter by Species
    if (is.null(input$chain_filter_species)){
      data <- data
    }  else if (any(input$chain_filter_species != "")){
      data = data[which(data$Species%in%input$chain_filter_species),]
    }
    return(data)
  })
  # For interaction_table
  curFilteredData_interaction_table <- reactive({
    # The full interaction_table is huge. So, by default we only load interactions related to user-selected chains.
    # So, depending on user selection, interaction_table is loaded directly from the MySQL database (using function interaction_table_by_chainid)
    #chain_ids = which(selection_all_chain$all$sel==T)
    #cat(chain_ids)
    #data <- ifelse (length(chain_ids)==0,
    #                interaction_table_show_all(), # If none chain is selected, select the complete table (slow)
    #                interaction_table_by_chainid(chain_ids[!is.na(chain_ids)]))
    data = interaction_table
    # Filter by PDB
    if (is.null(input$interaction_filter_pdb)){
      data <- data
    }  else if (any(input$interaction_filter_pdb != "")){
      data <- data[which(data$Pdb%in%input$interaction_filter_pdb),]
    }
    # Filter by Species
    if (is.null(input$interaction_filter_species)){
      data <- data
    }  else if (any(input$interaction_filter_species != "")){
      data = data[which(data$Species%in%input$interaction_filter_species),]
    }
  })
  #### ####
  
  #### Reset filtering button ####
  # complex_table
  observeEvent(input$reset_button, {
    shinyjs::reset("filter_pdb")
    shinyjs::reset("filter_experiment")
    shinyjs::reset("filter_resolution")
    shinyjs::reset("filter_species")
  })
  # chain_table
  observeEvent(input$reset_button_chain, {
    shinyjs::reset("chain_filter_pdb")
    shinyjs::reset("chain_filter_species")
  })
  # interaction_table
  observeEvent(input$reset_button_interaction, {
    shinyjs::reset("interaction_filter_pdb")
    shinyjs::reset("interaction_filter_species")
  })
  #### ####
  
  #################################################################################
  #### TABLE: Complex #############################################################
  #################################################################################
  
  #### Render the table containing shiny inputs ####
  output$tbl = DT::renderDataTable({
    # Get current data 
    data <- curFilteredData_complex_table()
    
    # Switch to show only selected entries 
    if (input$complex_selection_switch){
      data <- data[which(selection_all$all$sel==T),]
    }else{
      data <- data
    }
      
    # Selected entries
    selected_ids =  which(rownames(data)%in%which(selection_all$all$sel==T))
    
    # Convert by pretty colnames
    colnames(data)=complex_tbl_dic$desc
    # Columns to show
    by_complex = as.character(unlist(input$complex_tbl_columnSelection_by_complex))
    by_chain = as.character(unlist(input$complex_tbl_columnSelection_by_chain))
    by_ss = as.character(unlist(input$complex_tbl_columnSelection_by_ss))
    data = data[,c(by_complex,by_chain,by_ss)]
    
    # The table
    datatable(data, 
              style = 'bootstrap',
              rownames= FALSE,
              selection = list(mode = "multiple", target= 'row',
                               selected = selected_ids),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                deferRender = TRUE
              )
    )
  }, server = T)
  #### ####
  
  #### Downloadable csv of selected dataset ####
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ComplexesTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(curFilteredData_complex_table(), file, row.names = FALSE)
    }
  )
  #### ####
  
  #### Show Chains Button ####
  observeEvent(input$showChainsForSelection, {
    updateTabsetPanel(session, "menu", "chain")
  })
  #### ####
  
  #### Show Interactions Button ####
  observeEvent(input$showInteractionsForSelection, {
    updateTabsetPanel(session, "menu", "interaction")
  })
  #### ####
  
  #### Show Plot Buttons ####
  observeEvent(input$select_antibodyresidue_by_position_barplot, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_antibodyresidue_by_position_barplot')
  })
  observeEvent(input$select_atominteraction_residues_heatmap, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_atominteraction_residues_heatmap')
  })
  observeEvent(input$select_atominteraction_residues_network, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_network')
  })
  #### ####
  
  #################################################################################
  #### TABLE: Chain ###############################################################
  #################################################################################
  
  #### Render the table containing shiny inputs ####
  output$chain_tbl = DT::renderDataTable({
    # Get current data 
    data <- curFilteredData_chain_table()
    
    # Switch to show only selected entries 
    if (input$chain_selection_switch){
      data <- data[which(selection_all_chain$all$sel==T),]
    }else{
      data <- data
    }
    
    # Selected entries
    selected_ids =  which(rownames(data)%in%which(selection_all_chain$all$sel==T))
    
    # Convert to pretty colnames
    colnames(data)=chain_tbl_dic$desc
    # Columns to show
    by_chain = as.character(unlist(input$chain_tbl_columnSelection_by_chain))
    by_ss = as.character(unlist(input$chain_tbl_columnSelection_by_ss))
    data = data[,c(by_chain,by_ss)]
    
    # The table
    datatable(data, 
              style = 'bootstrap',
              rownames= FALSE,
              selection = list(mode = "multiple", target= 'row', 
                               selected = selected_ids),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                deferRender = TRUE
              )
    )
  }, server = T)
  #### ####
  
  #### Downloadable csv of selected dataset ####
  output$downloadData_chain <- downloadHandler(
    filename = function() {
      paste("ChainsTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(curFilteredData_chain_table(), file, row.names = FALSE)
    }
  )
  #### ####
  
  #### Show Complexes Button ####
  observeEvent(input$chain_showComplexesForSelection, {
    updateTabsetPanel(session, "menu", "complex")
  })
  #### ####
  
  #### Show Interactions Button ####
  observeEvent(input$chain_showInteractionsForSelection, {
    updateTabsetPanel(session, "menu", "interaction")
  })
  #### ####
  
  #### Show Plot Buttons ####
  observeEvent(input$chain_select_antibodyresidue_by_position_barplot, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_antibodyresidue_by_position_barplot')
  })
  observeEvent(input$chain_select_atominteraction_residues_heatmap, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_atominteraction_residues_heatmap')
  })
  observeEvent(input$chain_select_atominteraction_residues_network, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_network')
  })
  #### ####
  
  #################################################################################
  #### TABLE: Interaction #########################################################
  #################################################################################
  
  #### Render the table containing shiny inputs ####
  output$interaction_tbl = DT::renderDataTable({
    
    # Get current data 
    # If some chain was selected, shows only related interations
    chain_rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
    if (length(chain_rows_selected)==0){ # No chain selected
      data <- interaction_table_by_chainid(0)
    }else{ # Some chain selected
      data <- interaction_table_by_chainid(chain_rows_selected)
    }
    
    # Convert to pretty colnames
    colnames(data)=interaction_tbl_dic$desc
    # Columns to show
    by_antibody = as.character(unlist(input$interaction_tbl_columnSelection_by_antibody))
    by_antigen = as.character(unlist(input$interaction_tbl_columnSelection_by_antigen))
    by_interaction = as.character(unlist(input$interaction_tbl_columnSelection_by_interaction))
    data = data[,c(by_antibody,by_antigen,by_interaction)]
    
    # The table
    datatable(data, 
              style = 'bootstrap',
              rownames= FALSE,
              selection = list(mode = "multiple", target= 'row'),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                deferRender = TRUE
              )
    )
  }, server = T)
  #### ####
  
  #### Downloadable csv of selected dataset ####
  output$downloadData_interaction <- downloadHandler(
    filename = function() {
      paste("InteractionsTable", ".csv", sep = "")
    },
    content = function(file) {
      chain_rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
      if (length(chain_rows_selected)==0){ # No chain selected
        data <- interaction_table_by_chainid(0)
      }else{ # Some chain selected
        data <- interaction_table_by_chainid(chain_rows_selected)
      }
      write.csv(data, file, row.names = FALSE)
    }
  )
  #### ####
  
  #### Show Complexes Button ####
  observeEvent(input$interaction_showComplexesForSelection, {
    updateTabsetPanel(session, "menu", "complex")
  })
  #### ####
  
  #### Show Chains Button ####
  observeEvent(input$interaction_showChainsForSelection, {
    updateTabsetPanel(session, "menu", "chain")
  })
  #### ####
  
  #### Show Plot Buttons ####
  observeEvent(input$interaction_select_antibodyresidue_by_position_barplot, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_antibodyresidue_by_position_barplot')
  })
  observeEvent(input$interaction_select_atominteraction_residues_heatmap, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_atominteraction_residues_heatmap')
  })
  observeEvent(input$interaction_select_atominteraction_residues_network, {
    updateTabsetPanel(session, "menu", "interaction_charts")
    updateSelectInput(session, "int_plot_type",selected = 'select_network')
  })
  #### ####
  
  #################################################################################
  #### Plots: Stats ###############################################################
  #################################################################################
  
  #### Basic stats: Complexes ####
  number_of_complexes = nrow(unique(complex_table[,c("Pdb","HeavyChainName", "LightChainName", "AntigenChainName")]))
  number_of_pdbs = length(unique(complex_table[,"Pdb"]))
  number_of_antigenchains = sum(!is.na(unique(complex_table[,c("Pdb","AntigenChainName")])$AntigenChainName))
  number_of_heavychains = sum(!is.na(unique(complex_table[,c("Pdb","HeavyChainName")])$HeavyChainName))
  number_of_lightchains = sum(!is.na(unique(complex_table[,c("Pdb","LightChainName")])$LightChainName))
  dt = data.frame(
    row.names = c("Number of PBD structures:",
                  "Number of complexes (antibody+antigen):",
                  "Number of antigen chains:",
                  "Number of antibody heavy chains:",
                  "Number of antibody light chains:",
                  "Last release:"), 
    val = c(number_of_pdbs,
            number_of_complexes,
            number_of_antigenchains,
            number_of_heavychains,
            number_of_lightchains,
            "2019-03-25"),
    stringsAsFactors = F
  )
  output$general_stats_tbl = renderTable({dt},
                                         striped = TRUE,  
                                         spacing = 'xs',  
                                         width = '100%', align = 'l',
                                         rownames = TRUE, colnames = FALSE)
  #### ####
  
  #### Pdb by Experiment (Piechart) ####
  pdb_experiment = unique(complex_table[,c("Pdb","PDBExperiment")])
  output$pdbexperiment_piechart <- renderBillboarder({
    billboarder() %>% 
      bb_piechart(data = table(pdb_experiment[,"PDBExperiment"])) %>% 
      bb_labs(title = paste0("Number of PDBs by Experiment"))
  })
  #### ####
  
  #### Pdb Deposition - by Year (Barplot) ####
  pdb_deposition = unique(complex_table[,c("Pdb","PDBDepositionDate")])
  pdb_deposition$year = format(as.Date(pdb_deposition$PDBDepositionDate, format="%Y-%m-%d"),"%Y")
  pdb_deposition_table = table(pdb_deposition$year)
  pdb_deposition = as.data.frame(t(rbind(pdb_deposition_table,cumsum(pdb_deposition_table))))
  colnames(pdb_deposition) = c("deposited","accumulated")
  pdb_deposition$year = rownames(pdb_deposition)
  output$pdb_byyear_barplot <- renderBillboarder({
    billboarder() %>%
      bb_barchart(data = pdb_deposition[,c("year","accumulated","deposited")],stacked = TRUE) %>% 
      bb_data(names = list(deposited = "Deposited", accumulated = "Accumulated")) %>% 
      bb_y_grid(show = TRUE) %>%
      bb_x_axis(categories = pdb_deposition$year, fit = F, rotate = 90) %>% 
      bb_legend(position = "inset", inset = list(anchor = "top-left")) %>%
      bb_labs(title = "Number of PDBs by Year", y = "Number of structures") 
  })
  #### ####
  
  #### Pdb Resolution (Histogram) ####
  pdb_resolution = unique(complex_table[,c("Pdb","PDBResolution")])
  output$pdb_resolution_histogram <- renderBillboarder({
    billboarder() %>%
      bb_histogram(data = pdb_resolution$PDBResolution, binwidth = 0.3) %>%
      bb_data(names = list(y="Resolution")) %>%
      bb_x_axis(tick = list(values = seq(min(pdb_resolution$PDBResolution),max(pdb_resolution$PDBResolution)) ,fit = T)) %>%
      bb_labs(title = "Distribution of PDBs by Resolution", y = "Number of structures", x = "Resolution")
  })
  #### ####
  
  #### Complex Chains ####
  complex_chain = unique(complex_table[,c("Pdb","HeavyChainName","HeavyChainType",
                                          "LightChainName","LightChainType",
                                          "AntigenChainName","AntigenChainType")])
  #### ####
  
  #### antibody-producing organism ####
  antibody_species = unique(complex_table[,c("Pdb","HeavySpecies","LightSpecies")])
  
  heavy_species = as.data.frame(antibody_species$HeavySpecies, stringsAsFactors = F)
  #remove empty and NA 
  h_species = table(heavy_species[!apply(is.na(heavy_species) | heavy_species == "", 1, all),])
  h_species = c(sort(h_species,decreasing = T)[1:10],
                Others=sum(sort(h_species,decreasing = T)[11:length(h_species)]))
  
  light_species = as.data.frame(antibody_species$LightSpecies, stringsAsFactors = F)
  #remove empty and NA 
  l_species = table(light_species[!apply(is.na(light_species) | light_species == "", 1, all),])
  l_species = c(sort(l_species,decreasing = T)[1:10],
                Others=sum(sort(l_species,decreasing = T)[11:length(l_species)]))
  
  output$heavy_species_piechart <- renderBillboarder({
    billboarder() %>% 
      bb_piechart(data = as.table(h_species)) %>% 
      bb_labs(title = paste0("Antibody-producing organism (Heavy Chain)"))
  })
  
  output$light_species_piechart <- renderBillboarder({
    billboarder() %>% 
      bb_piechart(data = as.table(l_species)) %>% 
      bb_labs(title = paste0("Antibody-producing organism (Light Chain)"))
  })
  #### ####
  
  #### antigen-producing organism ####
  antigen_species = unique(complex_table[,c("Pdb","AntigenSpecies")])
  antigen_species = as.data.frame(antigen_species$AntigenSpecies, stringsAsFactors = F)
  #remove empty and NA 
  a_species = table(antigen_species[!apply(is.na(antigen_species) | antigen_species == "", 1, all),])
  a_species = c(sort(a_species,decreasing = T)[1:10],
                Others=sum(sort(a_species,decreasing = T)[11:length(a_species)]))
  output$antigen_species_piechart <- renderBillboarder({
    billboarder() %>% 
      bb_piechart(data = as.table(a_species)) %>% 
      bb_labs(title = paste0("Antigen-producing organism"))
  })
  #### ####
  
  #### Interactions by Residues ####
  sql = "SELECT ResidueInteraction.*, Rag.StructuralProperty, Rab.StructuralProperty FROM ResidueInteraction 
  LEFT JOIN Residue AS Rag ON ResidueInteraction.AntigenResidueId=Rag.Id
  LEFT JOIN Residue AS Rab ON ResidueInteraction.AntibodyResidueId=Rab.Id
  WHERE ?id"
  query <- sqlInterpolate(pool, sql, id = 1)
  ResidueInteraction_table = dbGetQuery(pool, query)
  
  interactions_by_residue = as.data.frame(rowsum(ResidueInteraction_table$Quantity,ResidueInteraction_table$Type))
  colnames(interactions_by_residue) = "Quantity"
  interactions_by_residue <- data.frame(Residues = row.names(interactions_by_residue), Quantity = interactions_by_residue$Quantity)
  
  output$interaction_by_residue_barplot <- renderBillboarder({
    billboarder() %>%
      bb_barchart(data = interactions_by_residue, stacked = F) %>% 
      bb_legend(show = FALSE) %>%
      bb_labs(x="Amino acid residue",y = "Number of interactions",title = "Quantity of interactions by residue")
  })
  #### ####

  #### Number of atom interactions by type ####
  output$atominteraction_by_type_barplot <- renderBillboarder({
    sql <- "SELECT 
          AtomInteraction.Type AS AtomInteraction_Type
          FROM AtomInteraction
          WHERE ?id"
    query <- sqlInterpolate(pool, sql, id = 1)
    interaction_table = dbGetQuery(pool, query)
    interaction_table$AtomInteraction_Type = merge(interaction_table,atom_interaction_type,by.x="AtomInteraction_Type",by.y="code")$type
    billboarder() %>%
      bb_barchart(data = table(interaction_table$AtomInteraction_Type), stacked = F) %>%
      bb_legend(show = F) %>%
      bb_labs(x = "Interaction type", y = "Count", title = "Number of atom interactions by type")
  })
  #### ####
  
  #### Distribution of atom interaction distances by type ####
  output$atominteraction_distances_by_type_histogram <- renderBillboarder({
    sql <- "SELECT 
          AtomInteraction.Type AS AtomInteraction_Type,
          AtomInteraction.Distance AS AtomInteraction_Distance
          FROM AtomInteraction
          WHERE ?id"
    query <- sqlInterpolate(pool, sql, id = 1)
    interaction_table = dbGetQuery(pool, query)
    interaction_table$AtomInteraction_Type = merge(interaction_table,atom_interaction_type,by.x="AtomInteraction_Type",by.y="code")$type
    billboarder() %>%
      bb_histogram(data = interaction_table, x = "AtomInteraction_Distance", group = "AtomInteraction_Type",
                   stacked = F) %>%
      bb_labs(x = "Interaction distance", y = "Count", title = "Distribution of atom interaction distances by type")
  })
  #### ####
  
  #### Distribution of atom interaction distances by Antigen amino acid ####
  output$atominteraction_distances_by_antigenresidue_histogram <- renderBillboarder({
    sql <- "SELECT 
            AtomInteraction.Distance AS AtomInteraction_Distance,
            AntigenResidue.AminoAcid             AS  AntigenResidue_AminoAcid
            FROM AtomInteraction
            LEFT JOIN Residue AS AntigenResidue ON AtomInteraction.AntigenResidueId=AntigenResidue.Id
            WHERE ?id"
    query <- sqlInterpolate(pool, sql, id = 1)
    interaction_table = dbGetQuery(pool, query)    
    billboarder() %>%
      bb_histogram(data = interaction_table, x="AtomInteraction_Distance", group = "AntigenResidue_AminoAcid",
                   stacked = F) %>%
      bb_labs(x = "Interaction distance", y = "Count", title = "Distribution of atom interaction distances by Antigen amino acid")
  })
  #### ####
  
  #### Distribution of atom interaction distances by Antibody amino acid ####
  output$atominteraction_distances_by_antibodyresidue_histogram <- renderBillboarder({
    sql <- "SELECT 
            AtomInteraction.Distance AS AtomInteraction_Distance,
            AntibodyResidue.AminoAcid             AS  AntibodyResidue_AminoAcid
            FROM AtomInteraction
            LEFT JOIN Residue AS AntibodyResidue ON AtomInteraction.AntibodyResidueId=AntibodyResidue.Id
            WHERE ?id"
    query <- sqlInterpolate(pool, sql, id = 1)
    interaction_table = dbGetQuery(pool, query)    
    
    billboarder() %>%
      bb_histogram(data = interaction_table, x = "AtomInteraction_Distance", group = "AntibodyResidue_AminoAcid",
                   stacked = F) %>%
      bb_labs(x = "Interaction distance", y = "Count", title = "Distribution of atom interaction distances by Antibody amino acid")
  })
  #### ####
  
  #################################################################################
  #### Plots: Interaction #########################################################
  #################################################################################
  
  #### Show information about current selected chains ####
  output$selection_info_complexTableView <- renderText({
    # Check for chains selected
    rows_selected_complex = (which(selection_all$all$sel==T))
    rows_selected = (which(selection_all_chain$all$sel==T))
    if (length(rows_selected)==0){ 
      return(("No chains selected!"))
    }else{
      # Get stats from selection (Chain as reference)
      data = curSelectedData_chain_table()
      nchains = length(unique(data$Id))
      nheavy = length(unique(data[which(data$ChainType=='H'),"Id"]))
      nlight = length(unique(data[which(data$ChainType=='L'),"Id"]))
      nantigen = length(unique(data[which(data$ChainType=='A'),"Id"]))
      ncomplex = length(complex_table[which(
        (complex_table$HeavyChainId%in%data$Id | complex_table$LightChainId%in%data$Id) &
                    complex_table$AntigenChainId%in%data$Id),"Id"])
      ninterac = nrow(interaction_table_by_chainid(data$Id))
      HTML(paste(
        (paste('Chains selected:',nchains)),
        (paste('Antibody Heavy chains:',nheavy)),
        (paste('Antibody Light chains:',nlight)),
        (paste('Antigen chains:',nantigen)),
        (paste('Complete Complexes:',ncomplex)),
        (paste('# Interactions:',ninterac)),
        sep="<br/>"))
    }
  })
  output$selection_info_chainTableView <- renderText({
    # Check for chains selected
    rows_selected_complex = (which(selection_all$all$sel==T))
    rows_selected = (which(selection_all_chain$all$sel==T))
    if (length(rows_selected)==0){ 
      return(("No chains selected!"))
    }else{
      # Get stats from selection (Chain as reference)
      data = curSelectedData_chain_table()
      nchains = length(unique(data$Id))
      nheavy = length(unique(data[which(data$ChainType=='H'),"Id"]))
      nlight = length(unique(data[which(data$ChainType=='L'),"Id"]))
      nantigen = length(unique(data[which(data$ChainType=='A'),"Id"]))
      ncomplex = length(complex_table[which(
        (complex_table$HeavyChainId%in%data$Id | complex_table$LightChainId%in%data$Id) &
          complex_table$AntigenChainId%in%data$Id),"Id"])
      ninterac = nrow(interaction_table_by_chainid(data$Id))
      HTML(paste(
        (paste('Chains selected:',nchains)),
        (paste('Antibody Heavy chains:',nheavy)),
        (paste('Antibody Light chains:',nlight)),
        (paste('Antigen chains:',nantigen)),
        (paste('Complete Complexes:',ncomplex)),
        (paste('# Interactions:',ninterac)),
        sep="<br/>"))
    }
  })
  output$selection_info_interactionTableView <- renderText({
    # Check for chains selected
    rows_selected_complex = (which(selection_all$all$sel==T))
    rows_selected = (which(selection_all_chain$all$sel==T))
    if (length(rows_selected)==0){ 
      return(("No chains selected!"))
    }else{
      # Get stats from selection (Chain as reference)
      data = curSelectedData_chain_table()
      nchains = length(unique(data$Id))
      nheavy = length(unique(data[which(data$ChainType=='H'),"Id"]))
      nlight = length(unique(data[which(data$ChainType=='L'),"Id"]))
      nantigen = length(unique(data[which(data$ChainType=='A'),"Id"]))
      ncomplex = length(complex_table[which(
        (complex_table$HeavyChainId%in%data$Id | complex_table$LightChainId%in%data$Id) &
          complex_table$AntigenChainId%in%data$Id),"Id"])
      ninterac = nrow(interaction_table_by_chainid(data$Id))
      HTML(paste(
        (paste('Chains selected:',nchains)),
        (paste('Antibody Heavy chains:',nheavy)),
        (paste('Antibody Light chains:',nlight)),
        (paste('Antigen chains:',nantigen)),
        (paste('Complete Complexes:',ncomplex)),
        (paste('# Interactions:',ninterac)),
        sep="<br/>"))
    }
  })
  #### ####
  
  #### Show information about current selected chains ####
  output$selection_info <- renderText({
    # Check for chains selected
    rows_selected = isolate(which(selection_all_chain$all$sel==T))
    if (length(rows_selected)==0){ 
      return(("No chains selected!"))
    }else{
      # Get stats from selection (Chain as reference)
      data = curSelectedData_chain_table()
      nchains = length(unique(data$Id))
      nheavy = length(unique(data[which(data$ChainType=='H'),"Id"]))
      nlight = length(unique(data[which(data$ChainType=='L'),"Id"]))
      nantigen = length(unique(data[which(data$ChainType=='A'),"Id"]))
      ncomplex = length(complex_table[which(
        (complex_table$HeavyChainId%in%data$Id | complex_table$LightChainId%in%data$Id) &
          complex_table$AntigenChainId%in%data$Id),"Id"])
      ninterac = nrow(interaction_table_by_chainid(data$Id))
      HTML(paste(
        (paste('Chains selected:',nchains)),
        (paste('Antibody Heavy chains:',nheavy)),
        (paste('Antibody Light chains:',nlight)),
        (paste('Antigen chains:',nantigen)),
        (paste('Complete Complexes:',ncomplex)),
        (paste('# Interactions:',ninterac)),
        sep="<br/>"))
    }
  })
  #### ####
  
  #### Generic interface (renderUI) to show different plot types (Plotly or Network) ####
  output$interaction_plots <- renderUI({
      switch(input$int_plot_type,
             "select_antibodyresidue_by_position_barplot" = plotlyOutput(outputId = "antibodyresidue_by_position_barplot"),
             "select_atominteraction_residues_heatmap" = plotlyOutput("atominteraction_residues_heatmap"),
             "select_network" = visNetworkOutput("network")
      )  
  })
  #### ####
  
  #### Interface to show plot datatables depending on user input (similar as above) ####
  output$interaction_plots_data <- renderUI({
    switch(input$int_plot_type,
           "select_antibodyresidue_by_position_barplot" = dataTableOutput("antibodyresidue_by_position_barplot_data"),
           "select_atominteraction_residues_heatmap" = dataTableOutput("atominteraction_residues_heatmap_data"),
           "select_network" = dataTableOutput("network_data")
    )
  })
  #### ####
  
  #### Frequency of Antibody amino acids by IMGT position (Plot) ####
  output$antibodyresidue_by_position_barplot <- renderPlotly({
    
    # Data - Check switch for selected-only
    # Check for chains selected
    rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
    if (input$interactionPlots_dataSelection_switch){
      tbl = residue_table_by_chainid(rows_selected)
    }else{
      tbl = residue_table_by_chainid(list())
    }
    data = curSelectedData_chain_table()
    
    # Get residue color scheme
    colorScheme = switch(input$resColor_imgtplot,
                         "by Amino acids" = 'AminoAcid',
                         "by Class" = 'class')
    legendLabel = switch(input$resColor_imgtplot,
                         "by Amino acids" = 'AA',
                         "by Class" = 'AA Class')
    
    # Get IMGT - IMGT postions definitions - http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
    choicePos = sapply(input$alignment_position_checkbox, switch,
                       "All" = 1:129,
                       "CDR1" = 27:38,
                       "CDR2" = 56:65,
                       "CDR3" = 105:117, # (105:116) for germline V-GENEs / (105:117) for rearranged V-J-GENEs or V-D-J-GENEs and corresponding cDNAs
                       "FR1" = 1:26,
                       "FR2" = 39:55,
                       "FR3" = 66:104,
                       "FR4" = 118:129
    )
    choicePos = unique(c(unlist(choicePos)))
    # Filter data by IMGT range
    tbl = tbl[which(tbl$IMGTIndex%in%choicePos),]
    
    # Get option for show by percentage
    barposition = switch(input$show_by_percentage_switch,
                         "Proportion (%)" = "fill",
                         "Raw counts" = "stack")
    ylabel = switch(input$show_by_percentage_switch,
                    "Proportion (%)" = "Frequency",
                    "Raw counts" = "Count")
    
    # Prepate table for ggplot
    tbl = plyr::count(tbl, c("IMGTIndex","AminoAcid"))
    tbl$IMGTIndex = as.factor(tbl$IMGTIndex)
    tbl = merge(tbl,aa_classes,by.x = "AminoAcid", by.y = "aa")
    
    # Plot
    p = ggplot(data=tbl, aes_string(x = "IMGTIndex", y = "freq", fill = colorScheme)) +
      geom_bar(position = barposition, stat = "identity") + 
      scale_fill_d3(palette = "category20") +
      labs(fill = legendLabel, x = "IMGT Index", y = ylabel) +
      theme_bw()
    ggplotly(p)
    
  })
  #### ####
  
  #### Frequency of Antibody amino acids by IMGT position (Data) ####
  output$antibodyresidue_by_position_barplot_data <- DT::renderDataTable({
    
    # Data - Check switch for selected-only
    # Check for chains selected
    rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
    if (input$interactionPlots_dataSelection_switch){
      tbl = residue_table_by_chainid(rows_selected)
    }else{
      tbl = residue_table_by_chainid(list())
    }
    
    # Get IMGT - IMGT postions definitions - http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
    choicePos = sapply(input$alignment_position_checkbox, switch,
                       "All" = 1:129,
                       "CDR1" = 27:38,
                       "CDR2" = 56:65,
                       "CDR3" = 105:117, # (105:116) for germline V-GENEs / (105:117) for rearranged V-J-GENEs or V-D-J-GENEs and corresponding cDNAs
                       "FR1" = 1:26,
                       "FR2" = 39:55,
                       "FR3" = 66:104,
                       "FR4" = 118:129
    )
    choicePos = unique(c(unlist(choicePos)))
    # Filter data by IMGT range
    tbl = tbl[which(tbl$IMGTIndex%in%choicePos),]
    
    # Prepate table for ggplot
    tbl = plyr::count(tbl, c("IMGTIndex","AminoAcid"))
    tbl$IMGTIndex = as.factor(tbl$IMGTIndex)
    tbl = merge(tbl,aa_classes,by.x = "AminoAcid", by.y = "aa")
    
    # The table
    datatable(tbl, 
              style = 'bootstrap',
              rownames= FALSE,
              extensions = 'Buttons', 
              selection = list(mode = "multiple", target= 'row'),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                deferRender = TRUE
              )
    )
  }, server = F)
  #### ####
  
  #### Number of interactions by antibody-antigen residues (Plot) ####
  output$atominteraction_residues_heatmap <- renderPlotly({
    
    # Data 
    # Check switch for selected-only
    if (input$interactionPlots_dataSelection_switch){
      # Check for chains selected
      rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
      if(length(rows_selected)>0){
        interaction_table = interaction_table_by_chainid(rows_selected)
      }else{
        interaction_table = interaction_table_show_all()
      }
    }else{
      interaction_table = interaction_table_show_all()
    }
    
    # Filter by type
    choiceType = input$interaction_type_picker
    if (length(choiceType)>0){
      df = interaction_table[which(interaction_table$AtomInteraction_Type%in%choiceType),]
    }else{
      df = data.frame(AntigenResidue_AminoAcid=0,AntibodyResidue_AminoAcid=0,AtomInteraction_Distance=0)
    }
    df = as.data.frame(table(df[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid")]))
    df = as.data.frame(acast(df, AntigenResidue_AminoAcid~AntibodyResidue_AminoAcid, value.var="Freq"))
    
    # Type of data
    if (input$interaction_heatmap_data == "Distance") {
      df = interaction_table[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid","AtomInteraction_Distance")]
      df = switch(input$interaction_heatmap_data_summary,
                  "Sum" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, sum)),
                  "Median" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, median)),
                  "Mean" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, mean))
      )
      df = xtabs(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df)
    }

    # Clustering and reordering rows and columns
    row.order <- hclust(dist(df))$order
    col.order <- hclust(dist(t(df)))$order
    df = df[row.order, col.order]    
    
    plot_ly(x=colnames(df),y=rownames(df),z=as.matrix(df), type = "heatmap",colorscale = "Greys") %>%
      colorbar(title = "Count") %>%
      layout(xaxis = list(title = "Antibody residues"),
             yaxis = list(title = "Antigen residues")
             #,title = "Number of interactions by antibody-antigen residues"
      )
  })
  #### ####
  
  #### Number of interactions by antibody-antigen residues (Data) ####
  output$atominteraction_residues_heatmap_data <- DT::renderDataTable({
    
    # Data 
    # Check switch for selected-only
    if (input$interactionPlots_dataSelection_switch){
      # Check for chains selected
      rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
      if(length(rows_selected)>0){
        interaction_table = interaction_table_by_chainid(rows_selected)
      }else{
        interaction_table = interaction_table_show_all()
      }
    }else{
      interaction_table = interaction_table_show_all()
    }
    
    # Filter by type
    choiceType = input$interaction_type_picker
    if (length(choiceType)>0){
      df = interaction_table[which(interaction_table$AtomInteraction_Type%in%choiceType),]
    }else{
      df = data.frame(AntigenResidue_AminoAcid=0,AntibodyResidue_AminoAcid=0,AtomInteraction_Distance=0)
    }
    df = as.data.frame(table(df[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid")]))
    df = as.data.frame(acast(df, AntigenResidue_AminoAcid~AntibodyResidue_AminoAcid, value.var="Freq"))
    
    # Type of data
    if (input$interaction_heatmap_data == "Distance") {
      df = interaction_table[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid","AtomInteraction_Distance")]
      df = switch(input$interaction_heatmap_data_summary,
                  "Sum" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, sum)),
                  "Median" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, median)),
                  "Mean" = as.data.frame(aggregate(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df, mean))
      )
      df = xtabs(AtomInteraction_Distance ~ AntigenResidue_AminoAcid+AntibodyResidue_AminoAcid, df)
    }
    
    # The table
    datatable(df, 
              style = 'bootstrap',
              rownames= T,
              extensions = 'Buttons', 
              selection = list(mode = "multiple", target= 'row'),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                deferRender = TRUE
              )
    )
    
  }, server = F)
  #### ####
  
  #### Network (Plot) ####
  output$network <- renderVisNetwork({
    
    # Data 
    # Check switch for selected-only
    if (input$interactionPlots_dataSelection_switch){
      # Check for chains selected
      rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
      if(length(rows_selected)>0){
        interaction_table = interaction_table_by_chainid(rows_selected)
      }else{
        interaction_table = interaction_table_by_chainid(23)
      }
    }else{
      interaction_table = interaction_table_by_chainid(23)
    }
    
    # Merge chain info with cath domains table
    cathtbl = cathdb[which(cathdb$PdbId%in%interaction_table$AntibodyChain_Pdb & cathdb$ChainId%in%interaction_table$AntibodyChain_ChainName),]
    
    # Filter by type
    choiceType = input$interacNet_type_picker
    if (length(choiceType)>0){
      interaction_table = interaction_table[which(interaction_table$AtomInteraction_Type%in%choiceType),]
    }else{
      interaction_table = data.frame(AntigenResidue_AminoAcid=0,AntibodyResidue_AminoAcid=0,AtomInteraction_Distance=0)
    }
    
    data = interaction_table[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid","AtomInteraction_Type")]
    typecolormap = as.data.frame(cbind(levels(data$AtomInteraction_Type),pal_d3("category20c")(18)))
    colnames(typecolormap) = c("Type","Color")
    data$color = typecolormap$Color[match(data$AtomInteraction_Type,typecolormap$Type)]
    data$AtomInteraction_Type = as.character(data$AtomInteraction_Type)
    
    freq_of_interaction = as.data.frame(table(data))
    freq_of_interaction = freq_of_interaction[which(!freq_of_interaction$Freq==0),]
    freq_of_interaction$AtomInteraction_Type = as.factor(freq_of_interaction$AtomInteraction_Type)
    
    #freq_of_interaction$Freq[which(freq_of_interaction$Freq<input$interacNet_slider)] = 0
    freq_of_interaction = freq_of_interaction[which(freq_of_interaction$Freq>=input$interacNet_slider),]
    
    node_ids = (c(
      paste0('Ag_',freq_of_interaction$AntigenResidue_AminoAcid),
      paste0('Ab_',freq_of_interaction$AntibodyResidue_AminoAcid)
    ))
    node_labels = (c(
      as.character(freq_of_interaction$AntigenResidue_AminoAcid),
      as.character(freq_of_interaction$AntibodyResidue_AminoAcid)
    ))
    node_colours = (c(
      rep("#D43F3AFF",length(freq_of_interaction$AntigenResidue_AminoAcid)),
      rep("#46B8DAFF",length(freq_of_interaction$AntibodyResidue_AminoAcid))
    ))
    node_group = (c(
      rep("Antigen",length(freq_of_interaction$AntigenResidue_AminoAcid)),
      rep("Antibody",length(freq_of_interaction$AntibodyResidue_AminoAcid))
    ))
    # removing redundancy
    unik = !duplicated(node_ids)  ## logical vector of unique values
    unik_indexes = seq_along(node_ids)[unik]  ## indices
    
    nodes <- data.frame(id = node_ids[unik_indexes],
                        label = node_labels[unik_indexes],
                        group = node_group[unik_indexes],
                        color = node_colours[unik_indexes],
                        shape='ellipse', value = 10,
                        title = paste0("<p>",aa_classes$class[match(node_labels[unik_indexes],aa_classes$aa)],"</p>")
                        )
    
    edges <- data.frame(from = paste0('Ag_',data$AntigenResidue_AminoAcid), 
                        to = paste0('Ab_',data$AntibodyResidue_AminoAcid),
                        color = data$color, width = 3,
                        title = paste0("<p>", data$AtomInteraction_Type,"</p>"))
    
    edges <- data.frame(from = paste0('Ag_',freq_of_interaction$AntigenResidue_AminoAcid), 
                        to = paste0('Ab_',freq_of_interaction$AntibodyResidue_AminoAcid),
                        color = freq_of_interaction$color, width = freq_of_interaction$Freq, 
                        label = freq_of_interaction$AtomInteraction_Type,
                        title = paste0("<p>", freq_of_interaction$AtomInteraction_Type,"<br>", freq_of_interaction$Freq, " interactions</p>"))

    # edges data.frame for legend
    ledges <- data.frame(color = typecolormap$Color)

    visNetwork(nodes, edges, width = "100%") %>%
      visGroups(groupname = "Antigen", color = "#D43F3AFF") %>%
      visGroups(groupname = "Antibody", color = "#46B8DAFF") %>%
      visNodes(scaling = list(label = list(enabled = T))) %>%
      visOptions(highlightNearest = list(enabled = T, degree = 2, hover = T), 
                 nodesIdSelection = TRUE,  
                 collapse = TRUE,
                 selectedBy = "group") %>% 
      visLayout(randomSeed = 123) %>%
      visLegend()
    
  })
  #### ####
  
  #### Network (Data) ####
  output$network_data <- DT::renderDataTable({
    
    # Data 
    # Check switch for selected-only
    if (input$interactionPlots_dataSelection_switch){
      # Check for chains selected
      rows_selected = selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)]
      if(length(rows_selected)>0){
        interaction_table = interaction_table_by_chainid(rows_selected)
      }else{
        interaction_table = interaction_table_by_chainid(23)
      }
    }else{
      interaction_table = interaction_table_by_chainid(23)
    }
    
    # Filter by type
    choiceType = input$interacNet_type_picker
    if (length(choiceType)>0){
      interaction_table = interaction_table[which(interaction_table$AtomInteraction_Type%in%choiceType),]
    }else{
      interaction_table = data.frame(AntigenResidue_AminoAcid=0,AntibodyResidue_AminoAcid=0,AtomInteraction_Distance=0)
    }
    
    data = interaction_table[,c("AntigenResidue_AminoAcid","AntibodyResidue_AminoAcid")]
    
    freq_of_interaction = as.data.frame(table(data))
    
    #freq_of_interaction$Freq[which(freq_of_interaction$Freq<input$interacNet_slider)] = 0
    freq_of_interaction = freq_of_interaction[which(freq_of_interaction$Freq>=input$interacNet_slider),]
    
    # The table
    datatable(freq_of_interaction, 
              style = 'bootstrap',
              rownames= FALSE,
              extensions = 'Buttons', 
              selection = list(mode = "multiple", target= 'row'),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                deferRender = TRUE
              )
    )
  }, server = F)
  #### ####
  
  #################################################################################
  #### Plots: Custom ##############################################################
  #################################################################################

  #### Columns to show as options ####
  complex_discrete_variables = c("AffinityTemperature")
  complex_continuous_variables = c("DeltaG","Affinity",
                                   "PDBResolution","PDBRFree",
                                   "PDBRFactor",
                                   "HeavyGravy","HeavyIsoelectricPoint","HeavyTotalNumberofAas",
                                   "HeavyAliphaticIndex","HeavyMolecularWeight",
                                   "HeavyNegativeCharged","HeavyPositiveCharged",             
                                   "HeavyPolarUncharged","HeavySmall",                        
                                   "HeavyHydrophobic","HeavyAlcohol",                   
                                   "HeavyAromatic","HeavyAlphaHelixStructure",   
                                   "HeavyIsolatedBetaBridgeStructure","HeavyStrandStructure",              
                                   "HeavyHelix3_10Structure","HeavyPiHelixStructure",            
                                   "HeavyTurnStructure","HeavyBendStructure",           
                                   "HeavyNoneStructure",
                                   "LightGravy","LightIsoelectricPoint","LightTotalNumberofAas",
                                   "LightAliphaticIndex","LightMolecularWeight",
                                   "LightNegativeCharged","LightPositiveCharged",             
                                   "LightPolarUncharged","LightSmall",                        
                                   "LightHydrophobic","LightAlcohol",                   
                                   "LightAromatic","LightAlphaHelixStructure",   
                                   "LightIsolatedBetaBridgeStructure","LightStrandStructure",              
                                   "LightHelix3_10Structure","LightPiHelixStructure",            
                                   "LightTurnStructure","LightBendStructure",           
                                   "LightNoneStructure",
                                   "AntigenGravy","AntigenIsoelectricPoint","AntigenTotalNumberofAas",
                                   "AntigenAliphaticIndex","AntigenMolecularWeight",
                                   "AntigenNegativeCharged","AntigenPositiveCharged",             
                                   "AntigenPolarUncharged","AntigenSmall",                        
                                   "AntigenHydrophobic","AntigenAlcohol",                   
                                   "AntigenAromatic","AntigenAlphaHelixStructure",   
                                   "AntigenIsolatedBetaBridgeStructure","AntigenStrandStructure",              
                                   "AntigenHelix3_10Structure","AntigenPiHelixStructure",            
                                   "AntigenTurnStructure","AntigenBendStructure",           
                                   "AntigenNoneStructure")
  
  chain_continuous_variables = c("Gravy","IsoelectricPoint","AliphaticIndex","MolecularWeight",
                                 "TotalNumberOfAas","NegativeCharged","PositiveCharged","PolarUncharged",
                                 "Small","Hydrophobic","Alcohol","Aromatic","AlphaHelixStructure",
                                 "IsolatedBetaBridgeStructure","StrandStructure","Helix3_10Structure",
                                 "PiHelixStructure","TurnStructure","BendStructure","NoneStructure")
  chain_discrete_variables = c("Pdb","ChainName","ChainType","Species")
  
  interaction_continuous_variables = c("AtomInteraction_Distance",
                                       "AntibodyResidue_SASComplexed","AntibodyResidue_SASUncomplexed",
                                       "AntigenResidue_SASComplexed","AntigenResidue_SASUncomplexed")
  interaction_discrete_variables = c("AtomInteraction_Type",
                                     "AtomInteraction_AntibodyAtom","AntibodyResidue_AminoAcid",
                                     "AntibodyResidue_StructuralProperty","AntibodyResidue_IMGTIndex",
                                     "AtomInteraction_AntigenAtom","AntigenResidue_AminoAcid",
                                     "AntigenResidue_StructuralProperty","AntigenResidue_IMGTIndex")
  #### ####
  
  #### Variable 1 ####
  output$variable <- renderUI({
    obj <- switch(input$dataset,
                "complex" = curSelectedData_complex_table(),
                "chain" = curSelectedData_chain_table(),
                "interaction" = interaction_table_by_chainid(selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)])
    )
    var.opts <- switch(input$dataset,
                     "complex" = complex_continuous_variables,
                     "chain" = chain_continuous_variables,
                     "interaction" = interaction_continuous_variables)
    str_label <- switch(input$plot_mode,
                       "One" = "Variable:",
                       "Two" = "Variable (y-axis):")
    selectInput("variable",str_label, var.opts) # uddate UI 				 
  })
  #### ####
  
  #### Variable 2 ####
  output$variable2 <- renderUI({ 
    obj <- switch(input$dataset,
                  "complex" = curSelectedData_complex_table(),
                  "chain" = curSelectedData_chain_table(),
                  "interaction" = interaction_table_by_chainid(selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)])
    )
    var.opts <- switch(input$dataset,
                     "complex" = complex_continuous_variables,
                     "chain" = chain_continuous_variables,
                     "interaction" = interaction_continuous_variables)
    selectInput("variable2","Variable (x-axis):", var.opts) # uddate UI 				 
  })
  #### ####
  
  #### Group ####
  output$group <- renderUI({ 
    var.opts <- switch(input$dataset,
                       "complex" = c("PDBExperiment",
                                     "HeavySpecies",
                                     "LightSpecies",
                                     "AntigenSpecies"),
                       "chain" = chain_discrete_variables,
                       "interaction" = interaction_discrete_variables)
    selectInput("group","Groups:",c(None='NA',var.opts),selectize=FALSE) # update UI 				 
  })
  #### ####
  
  #### Barplot options ####
  output$show_points <- renderUI({
    checkboxInput("show.points", "Show points", TRUE)
  })
  #### ####
  
  #### Plot type ####
  output$plot_type <- renderUI({
    if(input$plot_mode=="Two") {
      selectInput("plot.type","Plot Type:", 
                  list('Scatter Plot (points)' = "scatter_geom_bin", 
                       'Scatter Plot (hex bins)' = "scatter_geom_hex"
                  )
      )
    }else{
      selectInput("plot.type","Plot Type:", 
                  list('Histogram' = "histogram",
                       'Density Plot' = "density",
                       'Bar Plot' = "bar",
                       'Box Plot' = "boxplot"
                  )
      )
    }
  })
  #### ####
  
  #### Plot itself ####
  output$cplot <- renderPlotly({

    plot.obj <- list()
    plot.obj$data <- switch(input$dataset,
                            "complex" = curSelectedData_complex_table(),
                            "chain" = curSelectedData_chain_table(),
                            "interaction" = interaction_table_by_chainid(selection_all_chain$all$ID[which(selection_all_chain$all$sel==T)])
    )
    
    if (input$group=='NA'){
      plot.obj$group <- NA
    } else {
      plot.obj$group <- with(plot.obj$data,get(input$group))
    }
    
    plot.obj$variable <- with(plot.obj$data,get(input$variable))
    
    # Set plot type
    if(input$plot_mode == "One"){
      plot.type <- switch(input$plot.type,
                        "boxplot" 	= geom_boxplot(),
                        "histogram" =	geom_histogram(alpha=0.5, position="dodge"),
                        "density" 	=	geom_density(alpha=.75),
                        "bar" 		  =	geom_bar(position="dodge", stat="identity")
      )
    } else {
      plot.obj$variable2 <- with(plot.obj$data,get(input$variable2)) 
      #dynamic plotting options
      plot.type <- switch(input$plot.type,
                        "scatter_geom_bin" 	= geom_point(),
                        "scatter_geom_hex"  =	geom_hex()
      )
    }
    
    # By plot type...
    if (input$plot.type=="boxplot")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= input$variable, 
                       y 		= plot.obj$variable
                     )
        ) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$group, 
                       y 		= plot.obj$variable,
                       fill 	= as.factor(plot.obj$group)
                     )
        ) + plot.type
      }
      if(input$show.points==TRUE){ 
        p <- p + geom_point(color='black', alpha=0.5, position = 'jitter')
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= '',
        y 		= input$variable
      )  + theme_bw()
    }
    
    if (input$plot.type=="histogram")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, aes(x = plot.obj$variable)) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		  = plot.obj$variable,
                       fill 	= as.factor(plot.obj$group)
                     )
        ) + plot.type
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= input$variable,
        y 		= "Count"
      )  + theme_bw()
    }
    
    if (input$plot.type=="density")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, aes(x = plot.obj$variable)) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		  = plot.obj$variable,
                       fill 	= as.factor(plot.obj$group)
                     )
        ) + plot.type
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= input$variable,
        y 		= "Percent"
      )  + theme_bw()
    }
    
    if (input$plot.type=="bar")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= input$variable, 
                       y 		= plot.obj$variable
                     )
        ) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$group, 
                       y 		= plot.obj$variable,
                       fill 	= as.factor(plot.obj$group)
                     )
        ) + plot.type
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= '',
        y 		= input$variable
      )  + theme_bw()
    }
    
    if (input$plot.type=="scatter_geom_bin")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$variable2, 
                       y 		= plot.obj$variable
                     )
        ) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$variable2, 
                       y 		= plot.obj$variable,
                       colour 	= as.factor(plot.obj$group)
                     )
        ) + plot.type
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= input$variable2,
        y 		= input$variable
      )  + theme_bw()
    }
    
    if (input$plot.type=="scatter_geom_hex")	{ 
      if (is.na(plot.obj$group)) {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$variable2, 
                       y 		= plot.obj$variable
                     )
        ) + plot.type
      } else {
        p <- ggplot (plot.obj$data, 
                     aes(
                       x 		= plot.obj$variable2, 
                       y 		= plot.obj$variable
                     )
        ) + plot.type + facet_wrap(~plot.obj$group)
      }
      # Axis Labels
      p <- p + labs(
        fill 	= input$group,
        x 		= input$variable2,
        y 		= input$variable
      )  + theme_bw()
    }
    
    print(ggplotly(p))

  })
  #### ####
}
###################################################################################

###################################################################################
# User interface (Main function) ----
ui <- function(){

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
             )),
             fluidRow(box(
               title = "Filter by PDB IDs", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, collapsed = T,
               plotOutput("plot3", height = 250)
             )),
             fluidRow(box(
               title = "Filter by antibody-producing organism", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, collapsed = T,
               plotOutput("plot3", height = 250)
             )),
             fluidRow(box(
               title = "Filter by antigen-producing organism", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, collapsed = T,
               plotOutput("plot3", height = 250)
             )),
             fluidRow(box(
               title = "Filter by germline classification", status = "primary", solidHeader = TRUE,
               collapsible = TRUE, collapsed = T,
               plotOutput("plot3", height = 250)
             ))

      )
    )
  }
  
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
      
      #,box(
      #  title = "Interaction Statistics",
      #  status = "info",
      #  width = 12,
      #  collapsible = T, 
      #  collapsed = F,
      #  billboarderOutput(outputId = "interaction_by_residue_barplot"),
      #  billboarderOutput(outputId = "atominteraction_by_type_barplot"),
      #  billboarderOutput(outputId = "atominteraction_distances_by_type_histogram"),
      #  billboarderOutput(outputId = "atominteraction_distances_by_antigenresidue_histogram"),
      #  billboarderOutput(outputId = "atominteraction_distances_by_antibodyresidue_histogram")
      #)
      
    )
  }
  #### ####
  
  #### Interaction plots page ####
  interaction_charts = function(){
    
    fluidPage(
      
      #### opts ####
      theme = shinytheme("yeti"),
      useShinyjs(),
      #### ####
      
      fluidRow(
        
        #### Select type of plot ####
        column(12,                 
               selectInput("int_plot_type",h4("Plot by"), width = '100%',
                           list('Frequency of Antibody amino acids by IMGT position' = "select_antibodyresidue_by_position_barplot",
                                'Number of interactions by antibody-antigen residues' = "select_atominteraction_residues_heatmap",
                                'Network' = 'select_network')
               )
        )
        #### ####
        
      ),
      
      fluidRow(
        #### Plot area ####
        tabBox(
          title = "Explore interactions",
          id = "interactionBox",
          width = 8,
          height = 600,
          tabPanel("Plot",  withSpinner(uiOutput("interaction_plots"))),
          tabPanel("Data", withSpinner(uiOutput("interaction_plots_data")))
        ),
        #### ####
        
        #### General options for all kinds of charts ####
        box(
          title = "Selection info",
          status = "primary",
          collapsible = T,
          width = 4,
          materialSwitch(inputId = "interactionPlots_dataSelection_switch",
                        label = "Show only selected entries", 
                        status = "primary", right = F, value=T),
          withSpinner(htmlOutput("selection_info"))
        ),
        #### ####
        
        #### Plot options menu (Conditional) #### 
        box(
          title = "Plot options",
          status = "warning",
          width = 4,
          
          #### Options for plot by IMGT position ####
          conditionalPanel(condition = 'input.int_plot_type == "select_antibodyresidue_by_position_barplot"',
                           # Select IMGT position
                           pickerInput(
                             inputId = "alignment_position_checkbox", 
                             label = "Show by IMGT position:", 
                             choices = c("CDR1", "CDR2", "CDR3",
                                         "FR1","FR2","FR3","FR4"),
                             selected = "CDR1", 
                             options = list(
                               'dropupAuto' = FALSE,
                               'actions-box' = TRUE, 
                               size = 5
                             ), 
                             multiple = TRUE
                           ),
                           # Show by count | percentage
                           radioButtons(inputId = "show_by_percentage_switch",
                                       label = "Show frequency by",
                                       choices = c("Raw counts","Proportion (%)"),
                                       selected = "Raw counts"
                           ),
                           # Color by aa | class
                           radioButtons(inputId = "resColor_imgtplot",
                                       label = "Color residues",
                                       choices = c("by Class","by Amino acids"),
                                       selected = "by Class"
                           )
                           
          ),
          #### ####
          
          #### Options for interaction heatmap ####
          conditionalPanel(condition = 'input.int_plot_type == "select_atominteraction_residues_heatmap"',
                           # Select interaction type
                           pickerInput(
                             inputId = "interaction_type_picker", 
                             label = "Filter by Atom Interaction type:", 
                             choices = as.character(atom_interaction_type$type),
                             selected = as.character(atom_interaction_type$type), 
                             options = list(
                               'dropupAuto' = FALSE,
                               'actions-box' = TRUE, 
                               size = 5
                             ), 
                             multiple = TRUE
                           ),
                           selectInput(
                             inputId = "interaction_heatmap_data", 
                             label = "Show by:", 
                             choices = c("Number of Interactions",
                                         "Distance"),
                             selected = "Number of Interactions"
                           ),
                           conditionalPanel(condition = 'input.interaction_heatmap_data == "Distance"',
                                            radioButtons(
                                              inputId = "interaction_heatmap_data_summary", 
                                              label = "Summarized by:",
                                              choices = c("Sum","Median","Mean"),
                                              selected="Sum"
                                            )
                           )
          ),
          #### ####
          
          #### Options for interaction network ####
          conditionalPanel(condition = 'input.int_plot_type == "select_network"',
                           pickerInput(
                             inputId = "interacNet_type_picker", 
                             label = "Filter by Atom Interaction type:", 
                             choices = as.character(atom_interaction_type$type),
                             selected = as.character(atom_interaction_type$type), 
                             options = list(
                               'dropupAuto' = FALSE,
                               'actions-box' = TRUE, 
                               size = 5
                             ), 
                             multiple = TRUE
                           ),
                           sliderInput(
                             inputId = "interacNet_slider",
                             label = "Minimum number of interactions to show:",
                             min = 0, max = 50,
                             value = 0
                           )
          )
          #### ####
        )
      )
    )
    
  }
  #### ####
  
  #### Custom graphs page ####
  custom_charts = function(){
    
    fluidPage(
      
      #### opts ####
      theme = shinytheme("yeti"),
      useShinyjs(),
      #### ####
      
      sidebarLayout(position = 'right',
        sidebarPanel(width = 4, 
          
          h4("Plot options"),
          
          #### Data Filter/Selection ####
          selectInput("dataset","Data:", 
                      list(Complexes = "complex", Chains = "chain", Interactions = "interaction")),
          radioButtons("plot_mode", 
                       "Number of Continuous Variables",
                       choices = c("One", "Two"), 
                       inline = TRUE),
          #### ####
          
          uiOutput("variable"), 	# depends on dataset ( set by output$variable in  server.R)
          conditionalPanel(
            condition = "input.plot_mode == 'Two'",
            uiOutput("variable2")
          ),
          uiOutput("group"),  		# depends on dataset	( set by output$group in  server.R)
          uiOutput("plot_type"),
          radioButtons("dataset1", 
                       "Data to show:",
                       choices = c("All", "Selected"), 
                       inline = TRUE),
          
          # Options for boxplot (Conditional) 
          conditionalPanel(condition = 'input.plot_type == "Box Plot"',
                           uiOutput("show_points")
          )
        ),
        
        #### Plot area ####
        mainPanel(width = 8, tabBox(
          title = "Explore the data",
          id = "customplotBox",
          width = 12,
          tabPanel("Plot", withSpinner(plotlyOutput("cplot")))
        ))
        #### ####
      )

    )
  }
  #### ####
  
  #### Complex table page ####
  complex_table_view = function(){
    fluidPage(
      
      #### Opts ####
      theme = shinytheme("yeti"),
      useShinyjs(),
      #### ####
      
      #### Table ####
      fluidRow(
        column(12,
               
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
                        #,infoBox(title = div(style="font-size: small;font-weight: normal",htmlOutput("selection_info_complexTableView")),
                        #        div(style="font-size: small;font-weight: normal",materialSwitch(inputId = "complex_selection_switch",label = "Show only selected entries", status = "primary", right = F)),
                        #        icon = icon("info"), fill = T)
                 )
               )
               #### ####
               ,br(),
               #### Table and filter sidebar ####
               box(
                 title = "Complexes Explorer",
                 id = "complexExplorerBox",
                 collapsible = T,
                 width = 10,
                 withSpinner(dataTableOutput("tbl"))
               ),
               box(
                 title = "Filter data by",
                 id = "dataFilterBox_Complex",
                 collapsible = T,
                 width = 2,
                 selectizeInput("filter_pdb", "PDB Id:",
                                choices = as.character(PdbList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 selectizeInput("filter_experiment", "PDB Experiment:",
                                choices = as.character(ExperimentList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 selectizeInput("filter_species", "Species:",
                                choices = as.character(SpeciesList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 sliderInput("filter_resolution", "PDB Resolution", 
                             min = min(complex_table$PDBResolution), 
                             max = max(complex_table$PDBResolution), 
                             value = c(min(complex_table$PDBResolution), 
                                       max(complex_table$PDBResolution))),
                 actionButton("select_all_button", "Select all", width = "100%"),
                 actionButton("deselect_all_button", "Deselect all", width = "100%"),
                 actionButton("reset_button", "Reset", width = "100%")
               )
               #### ####  

        )
      )
    )
  }
  #### ####
  
  #### Chain table page ####
  chain_table_view = function(){
    fluidPage(
      
      #### Opts ####
      theme = shinytheme("yeti"),
      useShinyjs(),
      #### ####
      
      #### Table ####
      fluidRow(
        column(12,
               
               #### Row on top - (options) ####
               fixedRow(
                 column(12
                        ,div(style="display:inline-block",downloadButton("downloadData_chain", "Download",class = "btn btn-primary"))
                        ,div(style="display:inline-block",dropdownButton(inputId = "chain_columnsToShow", label = "Show columns", 
                                                                         circle = F, status = "primary",
                                                                         selectizeInput("chain_tbl_columnSelection_by_chain", label = h5("Chain"), multiple = T, 
                                                                                        choices = as.character(chain_tbl_dic[which(chain_tbl_dic$class=="Chain"),c("desc")]),
                                                                                        selected = c("PDB","Chain Name","Chain Type","Chain Species")),
                                                                         selectizeInput("chain_tbl_columnSelection_by_ss", label = h5("Secondary Structure"), multiple = T,
                                                                                        choices = as.character(chain_tbl_dic[which(chain_tbl_dic$class=="SS"),c("desc")]))
                        ))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block",actionButton("chain_showComplexesForSelection","See complexes",class = "btn btn-success"))
                        ,div(style="display:inline-block",actionButton("chain_showChainsForSelection","See chains",class = "btn btn-success"))
                        ,div(style="display:inline-block",actionButton("chain_showInteractionsForSelection","See interactions",class = "btn btn-success"))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block",dropdownButton(inputId = "chain_table_seegraphs", label = "Plot interactions",
                                                                         circle = F, status = "warning", icon = icon("bar-chart-o"), 
                                                                         actionButton("chain_select_antibodyresidue_by_position_barplot","Frequency of Antibody amino acids by IMGT position",class = "btn btn-warning",width = '100%'),
                                                                         actionButton("chain_select_atominteraction_residues_heatmap","Number of interactions by antibody-antigen residues",class = "btn btn-warning", width = '100%'),
                                                                         actionButton("chain_select_atominteraction_residues_network","Network",class = "btn btn-warning", width = '100%')))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block;text-align:right",dropdownButton(
                          tags$h4("Selection info"),
                          div(style="font-size: small;font-weight: normal",htmlOutput("selection_info_chainTableView")),hr(),
                          div(style="font-size: small;font-weight: normal",materialSwitch(inputId = "chain_selection_switch",label = "Show only selected entries", status = "primary", right = F)),
                          circle = TRUE, status = "info", icon = icon("info")
                        ))
                 )
               )
               #### ####
               ,br(),
               #### Table and filter sidebar ####
               box(
                 title = "Chains Explorer",
                 id = "chainExplorerBox",
                 collapsible = T,
                 width = 10,
                 withSpinner(dataTableOutput("chain_tbl"))
               ),
               box(
                 title = "Filter data by",
                 id = "dataFilterBox_Chain",
                 collapsible = T,
                 width = 2,
                 selectizeInput("chain_filter_pdb", "PDB Id:",
                                choices = as.character(PdbList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 selectizeInput("chain_filter_species", "Species:",
                                choices = as.character(SpeciesList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 actionButton("select_all_button_chain", "Select all", width = "100%"),
                 actionButton("deselect_all_button_chain", "Deselect all", width = "100%"),
                 actionButton("reset_button_chain", "Reset", width = "100%")
               )
               #### ####
               
        )
      )
    )
  }
  #### ####
  
  #### Interaction table page ####
  interaction_table_view = function(){
    fluidPage(
      
      #### Opts ####
      theme = shinytheme("yeti"),
      useShinyjs(),
      #### ####
      
      #### Table ####
      fluidRow(
        column(12,
               
               #### Row on top - (options) ####
               fixedRow(
                 column(12,div(style="display:inline-block",downloadButton("downloadData_interaction", "Download",class = "btn btn-primary"))
                        ,div(style="display:inline-block",dropdownButton(inputId = "interaction_columnsToShow", label = "Show columns", 
                                                                         circle = F, status = "primary",
                                                                         selectizeInput("interaction_tbl_columnSelection_by_antibody", label = h5("Antibody"), multiple = T, 
                                                                                        choices = as.character(interaction_tbl_dic[which(interaction_tbl_dic$class=="Antibody"),c("desc")]),
                                                                                        selected = c("Antibody Atom","Amino acid (Antibody)")),
                                                                         selectizeInput("interaction_tbl_columnSelection_by_antigen", label = h5("Antigen"), multiple = T, 
                                                                                        choices = as.character(interaction_tbl_dic[which(interaction_tbl_dic$class=="Antigen"),c("desc")]),
                                                                                        selected = c("Antigen Atom","Amino acid (Antigen)")),
                                                                         selectizeInput("interaction_tbl_columnSelection_by_interaction", label = h5("Interaction"), multiple = T, 
                                                                                        choices = as.character(interaction_tbl_dic[which(interaction_tbl_dic$class=="Interaction"),c("desc")]),
                                                                                        selected = c("Interaction Type","Interaction Distance"))
                        ))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block",actionButton("interaction_showComplexesForSelection","See complexes",class = "btn btn-success"))
                        ,div(style="display:inline-block",actionButton("interaction_showChainsForSelection","See chains",class = "btn btn-success"))
                        ,div(style="display:inline-block",actionButton("interaction_showInteractionsForSelection","See interactions",class = "btn btn-success"))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block",dropdownButton(inputId = "interaction_table_seegraphs", label = "Plot interactions",
                                                                         circle = F, status = "warning", icon = icon("bar-chart-o"), 
                                                                         actionButton("interaction_select_antibodyresidue_by_position_barplot","Frequency of Antibody amino acids by IMGT position",class = "btn btn-warning",width = '100%'),
                                                                         actionButton("interaction_select_atominteraction_residues_heatmap","Number of interactions by antibody-antigen residues",class = "btn btn-warning", width = '100%'),
                                                                         actionButton("interaction_select_atominteraction_residues_network","Network",class = "btn btn-warning", width = '100%')))
                        ,div(style="display:inline-block",HTML("<span>&#124;</span>"))
                        ,div(style="display:inline-block;text-align:right",dropdownButton(
                          tags$h4("Selection info"),
                          div(style="font-size: small;font-weight: normal",htmlOutput("selection_info_interactionTableView")),hr(),
                          div(style="font-size: small;font-weight: normal",materialSwitch(inputId = "interaction_selection_switch",label = "Show only selected entries", status = "primary", right = F)),
                          circle = TRUE, status = "info", icon = icon("info")
                        ))
                 )
               )
               #### ####
               ,br(),
               #### Table ####
               box(
                 title = "Interactions Explorer",
                 id = "interactionExplorerBox",
                 collapsible = T,
                 width = 10,
                 withSpinner(dataTableOutput("interaction_tbl"))
               ),
               box(
                 title = "Filter data by",
                 id = "dataFilterBox_Interaction",
                 collapsible = T,
                 width = 2,
                 selectizeInput("interaction_filter_pdb", "PDB Id:",
                                choices = as.character(PdbList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 selectizeInput("interaction_filter_species", "Species:",
                                choices = as.character(SpeciesList),
                                multiple = TRUE,
                                options = list(placeholder = 'All')),
                 actionButton("select_all_button_interaction", "Select all", width = "100%"),
                 actionButton("deselect_all_button_interaction", "Deselect all", width = "100%"),
                 actionButton("reset_button_interaction", "Reset", width = "100%")
               )
               #### ####
               
        )
      )
    )
  }
  #### ####
  
  
  #### Dashboard page layout parts ----
  # Include: Header + Sidebars + Body
  
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
                complex_table_view()
        ),
        tabItem(tabName = "chain",
                chain_table_view()
        ),
        tabItem(tabName = "interaction",
                interaction_table_view()
        ),
        tabItem(tabName = "custom_charts",
                custom_charts()
        ),  
        tabItem(tabName = "interaction_charts",
                interaction_charts()
        ), 
        tabItem(tabName = "dbstats_charts",
                dbstats_charts()
        )
      )
    )
  #### ####
  
  #### Page call (main function) ----
  dashboardPagePlus(title = "Ydb",
                    header, 
                    sidebar, 
                    body, 
                    collapse_sidebar = T)
}
###################################################################################

###################################################################################
# Run app ----
shinyApp(ui, server)
###################################################################################