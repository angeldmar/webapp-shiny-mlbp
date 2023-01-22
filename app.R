#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(reticulate)

use_virtualenv("./renv/python/virtualenvs/renv-python-3.10")

source_python("utils/SmilesFunctions.py")

con_modelo_ac = gzfile("./ac/modelo_ac.rds")
modelo_ac <- readRDS(con_modelo_ac)
con_modelo_ad = gzfile("./ad/modelo_ad.rds")
modelo_ad <- readRDS(con_modelo_ad)
con_modelo_ai = gzfile("./ai/modelo_ai.rds")
modelo_ai <- readRDS(con_modelo_ai)
con_modelo_am = gzfile("./am/modelo_am.rds")
modelo_am <- readRDS(con_modelo_am)
con_modelo_ao = gzfile("./ao/modelo_ao.rds")
modelo_ao <- readRDS(con_modelo_ao)

con_predictores_filtrados_ac = gzfile("./ac/predictores_filtrados_ac.rds")
predictores_filtrados_ac <- readRDS(con_predictores_filtrados_ac)
con_predictores_filtrados_ad = gzfile("./ad/predictores_filtrados_ad.rds")
predictores_filtrados_ad <- readRDS(con_predictores_filtrados_ad)
con_predictores_filtrados_ai = gzfile("./ai/predictores_filtrados_ai.rds")
predictores_filtrados_ai <- readRDS(con_predictores_filtrados_ai)
con_predictores_filtrados_am = gzfile("./am/predictores_filtrados_am.rds")
predictores_filtrados_am <- readRDS(con_predictores_filtrados_am)
con_predictores_filtrados_ao = gzfile("./ao/predictores_filtrados_ao.rds")
predictores_filtrados_ao <- readRDS(con_predictores_filtrados_ao)

con_trained_recipe_ac = gzfile("./ac/trained_recipe_ac.rds")
trained_recipe_ac <- readRDS(con_trained_recipe_ac)
con_trained_recipe_ad = gzfile("./ad/trained_recipe_ad.rds")
trained_recipe_ad <- readRDS(con_trained_recipe_ad)
con_trained_recipe_ai = gzfile("./ai/trained_recipe_ai.rds")
trained_recipe_ai <- readRDS(con_trained_recipe_ai)
con_trained_recipe_am = gzfile("./am/trained_recipe_am.rds")
trained_recipe_am <- readRDS(con_trained_recipe_am)
con_trained_recipe_ao = gzfile("./ao/trained_recipe_ao.rds")
trained_recipe_ao <- readRDS(con_trained_recipe_ao)


ui <- navbarPage(

  theme = shinytheme("slate"),
  "Nombre de pagina",
  tabPanel("Predicciones de bioactividades",
    fluidPage(
      fluidRow(
        sidebarLayout(
          sidebarPanel(
            textInput("user_smiles", "Introduce tu codigo Smiles"),
            actionButton("prediction_button", "Calcular predicciones", class = "btn-block btn-lg")
          ),
          mainPanel(
            tableOutput("prediction_output")
          )
        )
      ),
      fluidRow(
        column(2,
          ),
        column(8,
          imageOutput("smiles_image")
        ),
        column(2,
          )
      ),
    )
  ),
  navbarMenu("Acerca de la pagina",
    tabPanel("Funcionamiento", "Bla bla"),
    tabPanel("Documentos", "Articulo cientifico sin publicar, pagina de tesis"),
    tabPanel("Precision de modelos", "lorem ipsum")
  ),
    tabPanel("Acerca de los autores", "loren ipsum")
)
  
  


server <- function(input, output, session) {
  
  # hay que agregar algo al render  
  # output$prediction_output <- renderTable()
   
   
  
  generate_image <- eventReactive(input$prediction_button, {
    drawing_image(input$user_smiles)
       
    })
   
   output$smiles_image <- renderImage({
     
     generate_image()
     width  <- session$clientData$output_smiles_image_width
     list(
       #ver si funciona en python enviarlo a al folder img
       src = "./img/2D_smiles.png",
       contentType = "image/png",
       alt = "Imagen de molecula",
       width = width
     )
     
   }, 
   deleteFile = F
   )
   
}

# Run the application 
shinyApp(ui = ui, server = server)
