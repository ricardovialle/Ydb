library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  #################################################################################
  #### HOME: Selection Menus ######################################################
  #################################################################################
  
  output$all_pdb_checkbox_sel <- renderPrint({ input$all_pdb_checkbox })
  
})
