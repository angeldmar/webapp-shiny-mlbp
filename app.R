library(shiny)
library(shinythemes)
library(reticulate)
library(caret)
library(recipes)
library(RSNNS)
library(RRF)

use_virtualenv("./renv/python/virtualenvs/renv-python-3.11")

smilesf <- import("utils.SmilesFunctions")

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
# El objeto RDS en ad guardo el algoritmo en vez de las variables seleccionadas
predictores_filtrados_ad <- predictores_filtrados_ad$optVariables
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


bioactivity_prediction <- function(descriptors_dataframe, filter_values, preprocess_recipe, prediction_model, prediction_type, activity_name) {
  
  descriptors_dataframe <- bake(preprocess_recipe, new_data = descriptors_dataframe)
  descriptors_dataframe <- descriptors_dataframe[, filter_values] 
  descriptors_dataframe <- bind_rows(descriptors_dataframe)
  
  column_names <- c("Bioactividad", "% Verdadero", "% Falso", "Predicción")
  
  if (prediction_type == "prob") {
    prediction <- predict(prediction_model, newdata = descriptors_dataframe, type = prediction_type) %>% 
      mutate('class'=names(.)[apply(., 1, which.max)])
    prediction_dataframe <- data.frame(bioactividad = activity_name, prediction)
    colnames(prediction_dataframe) <- column_names
    prediction_dataframe[c("% Verdadero", "% Falso")] <- lapply(prediction_dataframe[c("% Verdadero", "% Falso")], function(x) x * 100)
    return(prediction_dataframe)
  } else if (prediction_type == "raw") {
    prediction <- predict(prediction_model, newdata = descriptors_dataframe, type = prediction_type)
    prediction_dataframe <- data.frame(activity_name, "---", "---", as.character(prediction))
    colnames(prediction_dataframe) <- column_names
    return(prediction_dataframe)
  } else {
    return(print("Debe ingresar un tipo de predicción válida"))
  }
}


generate_predictions <- function(descriptors) {
  
  ac_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ac, trained_recipe_ac, modelo_ac, "prob", "Anticancerígeno")
  ad_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ad, trained_recipe_ad, modelo_ad, "raw", "Antidiabético")
  ai_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ai, trained_recipe_ai, modelo_ai, "prob", "Antiinflamatorio")
  am_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_am, trained_recipe_am, modelo_am, "raw", "Antimicrobiano")
  ao_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ao, trained_recipe_ao, modelo_ao, "prob", "Antioxidante")
  
  unified_predictions <- merge(x = ac_prediction, y = ad_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = ai_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = am_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = ao_prediction, all = TRUE)
  
  unified_predictions[unified_predictions == FALSE] <- "No posee bioactividad"
  unified_predictions[unified_predictions == TRUE] <- "Si posee bioactividad"
  
  return (unified_predictions)
}


ui <- navbarPage(

  theme = shinytheme("darkly"),
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
            tableOutput("prediction_table")
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

  user_predictions <- eventReactive(input$prediction_button, {
    user_descriptors_dataframe <- smilesf$generating_descriptors(input$user_smiles)
    predictions <- generate_predictions(user_descriptors_dataframe)
  })
  
  output$prediction_table <- renderTable({
    user_predictions()
    },     
    striped = TRUE,
    hover = TRUE,
    spacing = 'l',
    width = '100%',
    digits = 2,
    na = 'missing',
    align = 'c'
    )
  
  generate_image <- eventReactive(input$prediction_button, {
    smilesf$drawing_smiles(input$user_smiles)
    })
   
  output$smiles_image <- renderImage({
      
      generate_image()
      # width  <- session$clientData$output_smiles_image_width
      list(
        src = "./img/2D_smiles.png",
        contentType = "image/png",
        alt = "Imagen de molecula",
        width = '100%'
      )     
    }, 
    deleteFile = F
  ) 
}


# Run the application 
shinyApp(ui = ui, server = server)
