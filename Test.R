library(dplyr)
library(reticulate)
#library(RSNNS)
library(recipes)
library(caret)
library(RRF)
library(data.table) # solo si es necesario
library(randomForest)
library(inTrees)
library(snn)
library(wsrf)
library(stringr)

use_virtualenv("./renv/python/virtualenvs/renv-python-3.11")

# source_python("utils/SmilesFunctions.py")

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
    prediction_dataframe[c("% Verdadero", "% Falso")] <- lapply(prediction_dataframe[c("% Verdadero", "% Falso")], function(x) round(x * 100, digits = 2))
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

molecula_test<- "CN=C=O"

descriptors <- smilesf$generating_descriptors(molecula_test)
predictions_dataframe <- generate_predictions(descriptors)
# ver si es necesario
predictions_table <- as.data.table(predictions_dataframe)

# para shiny es
# renderTable(predictions_table)



# generating_descriptors <- function(molecule) {
#   # Hay que hacer una funcion de cleaning_input
#   descriptors2d <- smilesf$generating_descriptors2d_dataframe(molecule)
#   descriptors3d <- smilesf$generating_descriptors3d_dataframe(molecule)
#   lipinski_descriptors <- smilesf$calculating_lipinski_descriptors(molecule)
#   descriptors_dataframe <- merge(descriptors2d, 
#                                  descriptors3d, 
#                                  by.x = 0, 
#                                  by.y = 0,
#                                  all.x = TRUE,
#                                  all.y = TRUE)
#   descriptors_dataframe <- merge(descriptors_dataframe, 
#                                  lipinski_descriptors, 
#                                  by.x = "Row.names", 
#                                  by.y = 0,
#                                  all.x = TRUE,
#                                  all.y = TRUE)
#   rownames(descriptors_dataframe) <- descriptors_dataframe$Row.names
#   descriptors_dataframe <- descriptors_dataframe[!grepl("Row.names", names(descriptors_dataframe))]
#   
#   return(descriptors_dataframe)
# }

str_detect("casa", "c")
str_detect("gato", "c")
input_smiles <- "CN=C=O"
str_detect(input_smiles, "/^([^J][0-9BCOHNSOPrIFla@+\-\[\]\(\)\\\/%=#$]{6,})$/ig")
