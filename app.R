library(shiny)
library(shinythemes)
library(shinyFeedback)
library(reticulate)
library(stringr)
library(caret)
library(recipes)
library(RRF)
library(RSNNS)
library(snn)
library(inTrees)
library(wsrf)
library(randomForest)

# use_virtualenv("./renv/python/virtualenvs/renv-python-3.11")

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

bioactivity_prediction <- function(descriptors_dataframe, filter_values, preprocess_recipe, prediction_model, activity_name) {
  
  descriptors_dataframe <- bake(preprocess_recipe, new_data = descriptors_dataframe)
  descriptors_dataframe <- descriptors_dataframe[, filter_values] 
  descriptors_dataframe <- bind_rows(descriptors_dataframe)
  
  column_names <- c("Bioactividad"," Predicción")
  prediction <- predict(prediction_model, newdata = descriptors_dataframe, type = "raw")
  prediction_dataframe <- data.frame(activity_name, as.character(prediction))
  colnames(prediction_dataframe) <- column_names
  return(prediction_dataframe)
}


generate_predictions <- function(descriptors, id) {
  
  notify("Generando predicciones ...", id = id)
  on.exit(removeNotification(id), add = TRUE)
  
  ac_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ac, trained_recipe_ac, modelo_ac, "Anticancerígeno")
  ad_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ad, trained_recipe_ad, modelo_ad, "Antidiabético")
  ai_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ai, trained_recipe_ai, modelo_ai, "Antiinflamatorio")
  am_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_am, trained_recipe_am, modelo_am, "Antimicrobiano")
  ao_prediction <- bioactivity_prediction(descriptors, predictores_filtrados_ao, trained_recipe_ao, modelo_ao, "Antioxidante")
  
  unified_predictions <- merge(x = ac_prediction, y = ad_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = ai_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = am_prediction, all = TRUE)
  unified_predictions <- merge(x = unified_predictions, y = ao_prediction, all = TRUE)
  
  unified_predictions[unified_predictions == FALSE] <- "No"
  unified_predictions[unified_predictions == TRUE] <- "Sí"
  
  return (unified_predictions)
}


notify <- function(msg, id = NULL){
  showNotification(msg, id = id, duration = NULL, closeButton = FALSE)
}


ui <- navbarPage(
  
  theme = shinytheme("united"),
  "MLBioPrediction",
  tabPanel("Predicciones de bioactividades",
    fluidPage(
      tags$head(
        tags$style(HTML("
          
          a {
          display: block;
          }
          
          img {
            display: block;
            margin-left: auto;
            margin-right: auto;
            aspect-ratio: 16 / 9;
            max-height: 70vh;
            max-width: 100%;
            width: auto;
            position: relative;
            overflow: hidden;
          }
          
          .shiny-image-output {
            display: block;
            max-width: 100%;
            max-height: 100%;
            position: absolute;
            top: 0;
            bottom: 0;
            left: 0;
            right: 0;
          }
            
          .warning-box {
            text-align: center; 
            border-radius: 5px; 
            padding: 10px; 
            margin-bottom: 10px;
            margin-left: 5%;
            margin-right: 5%;
            cursor: pointer;
            color: white;
            background-color: #e95420;
          }
          
          .warning-title {
            margin-top: 0px;
            text-align: center; 
          }
        /*Es para cambiar de color con el hover,
        .link-box:hover {
          background-color: #ac3911;          
        }*/
          ")),

        tags$link( 
          rel = "shortcut icon", 
          type = "image/png", 
          href = "logo.png" 
        ), 
      ),
      useShinyFeedback(),
      fluidRow(
        sidebarLayout(
          sidebarPanel(
            textInput("user_smiles", "Introduce tu código SMILES canónico"),
            actionButton("prediction_button", "Calcular predicciones", class = "btn-block btn-lg")
          ),
          mainPanel(
            conditionalPanel(
              condition = "input.prediction_button == 0",
              
              tags$h2(class = "warning-title", "Limitaciones"),
              
              tags$div(
                class = "warning-box",
                "Es posible que en ocasiones genere información incorrecta.",
              ),
              tags$div(
                class = "warning-box",
                "Es probable que proporcione resultados más precisos al trabajar 
                con moléculas de origen natural y sus derivados, en comparación 
                con aquellas de origen sintético."
              ),
              tags$div(
                id = "box1",
                class = "warning-box link-box",
                "En la sección de Funcionamiento encontrarás información 
                detallada sobre los modelos y su precisión."
              ),
              tags$div(
                id = "box_2",
                class = "warning-box link-box",
                "Si necesitas información adicional y en mayor profundidad, 
                puedes consultar la documentación disponible en la página."
              ),
            ),
            conditionalPanel(
              condition = "input.prediction_button > 0",
              tableOutput("prediction_table")
            )
          ),
        )
        
      ),
      fluidRow(
        column(12,
          imageOutput("smiles_image")
        ),
      )
    )
  ),
  navbarMenu("Acerca de la página",
    tabPanel("Funcionamiento",
             id = "func",
             tags$h5("Esta página está diseñada para predecir actividades biológicas de estructuras 
                     moleculares utilizando el código SMILES canónico. Las predicciones se basan en 
                     modelos entrenados a partir de moleculas de origen natural y derivados a través 
                     de algoritmos de aprendizaje automatizado. "),
             tags$h5("A continuación, se describen los modelos utilizados y las precisiones obtenidas 
                     sobre moléculas de prueba, para la predicción de cada una de las bioactividades: "),
             tableOutput("model_info")
             ),
    tabPanel("Documentos",
             tags$h3("A continuación se presentan enlaces a documentos de posible interes:"),
             tags$div(
               HTML(# <a href="" target="_blank" rel="noopener noreferrer">Articulo cientifico sin publicar (opcional)</a>
                    # <a href="" target="_blank" rel="noopener noreferrer">Tesis de grado(si la subimos a algun lado)</a>
                    '<a href="https://github.com/angeldmar/Tesis-prediccion-bioactividades" target="_blank" rel="noopener noreferrer">Repositorio de los modelos de machine learning</a>
                    <a href="https://github.com/angeldmar/webapp-shiny-mlbp" target="_blank" rel="noopener noreferrer">Repositorio de la aplicación web</a>'
                  )
             )
    )
  ),
  tabPanel("Acerca de los autores", 
           tags$h5("Somos Angel Martínez, Erleigh Hogan y  Galilea Wug, Químicos Farmacéuticos egresados 
                   de la Facultad de Ciencias Químicas y Farmacia de la Universidad de San Carlos de Guatemala."), 
           tags$h5("Decidimos realizar este proyecto como investigadores del Departamento de Química Medicinal 
                    con la asesoría de nuestros catedráticos, MSc. Nereida Marroquín y Lic. Allan Vásquez, ya 
                    que reconocemos que la industria farmacéutica consume muchos recursos para encontrar moléculas 
                    con potencial biológico que puedan ser utilizadas para el desarrollo de nuevos fármacos.
                    "),
           tags$h5("Con esta herramienta web de libre acceso pretendemos ayudar en la reducción de recursos 
                   para que el investigador, solo utilizando el código SMILES, pueda predecir la actividad 
                   biológica de su estructura molecular para posteriormente confirmar dicha predicción a 
                   través de ensayos químicos específicos.")
  )
)
  
  
server <- function(input, output, session) {
  
  observeEvent(input$prediction_button, {
    updateTabsetPanel(session, "tabset", selected = "prediction_table")
  })
  
  clean_input <- eventReactive(input$prediction_button, {
    req(input$user_smiles)
    input_smiles <- gsub(" ", "", input$user_smiles)
    is_smile <-smilesf$validating_smiles(input_smiles)
    shinyFeedback::feedbackDanger("user_smiles", !is_smile, "Por favor introduzca un código SMILES válido")
    req(is_smile, cancelOutput = TRUE)
    return(input_smiles)
  })
  
  
  user_predictions <- reactive({
    id <- notify("Calculando descriptores...")
    on.exit(removeNotification(id), add = TRUE)
    user_descriptors_dataframe <- smilesf$generating_descriptors(clean_input())
    predictions <- generate_predictions(user_descriptors_dataframe, id)
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
  
  generate_image <- reactive({
    smilesf$drawing_smiles(clean_input())
    })
   
  output$smiles_image <- renderImage({
      
      generate_image()
      list(
        src = "./img/2D_smiles.png",
        contentType = "image/png",
        alt = "Imagen de molecula"
      )     
    }, 
    deleteFile = T
  ) 
  
  output$model_info <- renderTable({
    model_table_info = data.frame(x1 = c("Anticancerígena", "Antidiabética", 
                                         "Antiinflamatoria", "Antimicrobiano", 
                                         "Antioxidante"),  
                                  # x2 = c("Bosque aleatorio regularizado", 
                                  #        "Modelo basado en reglas de bosques aleatorios", 
                                  #        "Perceptrón multicapa", 
                                  #        "Clasificador de vecino más cercano estabilizado",
                                  #        "Bosque aleatorio de subespacio ponderado"),
                                  x2 = c("Regularized Random Forest", 
                                         "Random Forest Rule-Based Model", 
                                         "Multi-Layer Perceptron", 
                                         "Stabilized Nearest Neighbor Classifier",
                                         "Weighted Subspace Random Forest"),
                                  x3 = c("73.33%", "75.82%", "74.73%", "81.32%", "83.52%"))
    
    colnames(model_table_info) <- c("Bioactividad", "Modelo utilizado", "Precisión")
    return(model_table_info)
  },
    striped = TRUE,
    hover = TRUE,
    spacing = 'l',
    width = '100%',
    digits = 2,
    na = 'missing',
    align = 'c'
  )
}


# Run the application 
shinyApp(ui = ui, server = server)

# Para generar el reactlog
# reactlog::reactlog_enable() # antes de ejecutar la app
# shiny::reactlogShow() # despues de ejecutar la app


