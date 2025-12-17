# ====================================
# 0. Instalar/llamar paquetes
# ====================================
library(FinTS)
library(strucchange)
library(dynamac)
library(readxl)
library(urca)      # Para pruebas ADF y Johansen
library(tsDyn)     # Para VECM
library(dLagM)     # Para ARDL (o package "ARDL" en CRAN)
library(tseries)   # Pruebas ADF
library(dplyr)
library(ARDL)
library(ggplot2)
library(lmtest)
library(car)
library(bootUR)
library(stringr)
library(xtable)
# ====================================
# 1. Importar datos
# ====================================
df <- read_excel(file.choose(), sheet = "Hoja1")

# Variables esperadas:
# "Year", "ECI_SITC", "ECI_HS92", "FDI_GDP", "Trade_GDP", "GDP_capita",
# "GFCF_GDP", "School_enrollment", "GovQualityIndex", "ICT_Index"

df <- df %>%
  select(Year, ECI_SITC, ECI_HS92,FDI_GDP, Trade_GDP, GDP_capita, 
         GFCF_GDP, School_enrollment, GovQualityIndex, ICT_Index) %>%
  arrange(Year)

# Para convertir a objeto ts (time series anual)
# Como parte de replicar un paper, a menudo se usa ts() con start= c(1995)
df_ts <- ts(df[,-1], start=c(1995), frequency=1)

# ====================================
# 2. Exploración inicial
# ====================================
summary(df_ts)
# Plot de cada serie
plot(df_ts, main="Series Chile 1995-2023")
plot(df_ts[,1], main="ECI SITC 1995-2023", xlab="Year", ylab="ECI SITC")
plot(df_ts[,2], main="ECI HS92 1995-2023", xlab="Year", ylab="ECI HS92")
plot(df_ts[,3], main="Foreign Direct Investment (% of GDP) 1995-2023", xlab="Year", ylab="FDI as % of GDP")
plot(df_ts[,4], main="Trade (% of GDP) 1995-2023", xlab="Year", ylab="Trade (% of GDP)")
plot(df_ts[,5], main="GDP per Capita (current $) 1995-2023", xlab="Year", ylab="GDP per Capita (current $)")
plot(df_ts[,6], main="Gross Fixed Capital Formation (% of GDP) 1995-2023", xlab="Year", ylab="Gross Fixed Capital Formation (% of GDP)")
plot(df_ts[,7], main="Secondary School Enrollment % Gross 1995-2023", xlab="Year", ylab="Secondary School Enrollment % Gross")
plot(df_ts[,8], main="Government Quality Index 1995-2023", xlab="Year", ylab="Government Quality Index")
plot(df_ts[,9], main="Information and Communication Technology Index 1995-2023", xlab="Year", ylab="ICT Index")
# ====================================
# 3. Pruebas de estacionariedad
# ====================================

## ADF y PP para cada variable en I(0) / test adf con inteceptt o trend / PP con model constant o trend y fijo en z-tau




##I(0)
#Definir tablas para guardar resultados por variable y tipo de test y integracion
resultados_adf_I0 <- data.frame(
  Variable = character(),
  Tipo_Modelo = character(),
  Estadístico_ADF = character(),
  cval1 = numeric(),
  cval5 = numeric(),
  cval10 = numeric(),
  stringsAsFactors = FALSE
)

resultados_pp_I0 <- data.frame(
  Variable = character(),
  Tipo_Modelo = character(),
  Estadístico_PP = character(),
  cval1 = numeric(),
  cval5 = numeric(),
  cval10 = numeric(),
  stringsAsFactors = FALSE
)
#ADF
for (i in 1:ncol(df_ts)) {
  cat("\nVariable:", colnames(df_ts)[i])
  
  # datos variable
  x <- df_ts[, i]
  
  # ADF Test 
  adf_type = c("drift", "trend")
  
  for (model_type in adf_type){
    
    adf <- ur.df(x, type = model_type, selectlags = "BIC")
    adf_stat <- adf@teststat[1]
    if (model_type == 'drift') {
      rechazo_1 <- adf_stat < adf@cval[1]
      rechazo_5 <- adf_stat < adf@cval[3]
      rechazo_10 <- adf_stat < adf@cval[5]
      cval1 = adf@cval[1]
      cval5 = adf@cval[3]
      cval10 = adf@cval[5]
    }
    else {
      rechazo_1 <- adf_stat < adf@cval[1]
      rechazo_5 <- adf_stat < adf@cval[4]
      rechazo_10 <- adf_stat < adf@cval[7]
      cval1 = adf@cval[1]
      cval5 = adf@cval[4]
      cval10 = adf@cval[7]
    }
    
    if (rechazo_1) {
      resultado = sprintf("%f***", adf_stat)
    } 
    else if (rechazo_5) {
      resultado = sprintf("%f**", adf_stat)
    } 
    else if (rechazo_10) {
      resultado = sprintf("%f*", adf_stat)
    } 
    else {
      resultado = sprintf("%f", adf_stat)
    }
    
    #agregamos resultados
    resultados_adf_I0 <- rbind(resultados_adf_I0, data.frame(
      Variable = colnames(df_ts)[i],
      Tipo_Modelo = model_type,
      Estadístico_ADF = resultado,
      cval1 = cval1,
      cval5 = cval5,
      cval10 = cval10,
      stringsAsFactors = FALSE
    ))
    
  }
  
  pp_type = c("constant", "trend")
  # UR PP Test (constant y short)
  
  for (model_type in pp_type){
    pp <- ur.pp(x, type = "Z-tau", model = model_type, lags = 'short')
    stat_pp <- pp@teststat[1]
    
    crit_pp_1 <- pp@cval[1]
    crit_pp_5 <- pp@cval[2]
    crit_pp_10 <- pp@cval[3]
    
    
    rechazo_1 <- stat_pp < crit_pp_1
    rechazo_5 <- stat_pp < crit_pp_5
    rechazo_10 <- stat_pp < crit_pp_10
    if (rechazo_1) {
      resultado_pp = sprintf("%f***", stat_pp)
    } 
    else if (rechazo_5) {
      resultado_pp = sprintf("%f**", stat_pp)
    } 
    else if (rechazo_10) {
      resultado_pp = sprintf("%f*", stat_pp)
    } 
    else {
      resultado_pp = sprintf("%f", stat_pp)
    }
    
    #Agregar resultados a pp
    #agregamos resultados
    resultados_pp_I0 <- rbind(resultados_pp_I0, data.frame(
      Variable = colnames(df_ts)[i],
      Tipo_Modelo = model_type,
      Estadístico_PP = resultado_pp,
      cval1 = pp@cval[1],
      cval5 = pp@cval[2],
      cval10 = pp@cval[3],
      stringsAsFactors = FALSE
    ))
    
    
  }
}



##I(1)
#Definir tablas para guardar resultados por variable y tipo de test y integracion
resultados_adf_I1 <- data.frame(
  Variable = character(),
  Tipo_Modelo = character(),
  Estadístico_ADF = character(),
  cval1 = numeric(),
  cval5 = numeric(),
  cval10 = numeric(),
  stringsAsFactors = FALSE
)

resultados_pp_I1 <- data.frame(
  Variable = character(),
  Tipo_Modelo = character(),
  Estadístico_PP = character(),
  cval1 = numeric(),
  cval5 = numeric(),
  cval10 = numeric(),
  stringsAsFactors = FALSE
)
#ADF
for (i in 1:ncol(df_ts)) {
  cat("\nVariable:", colnames(df_ts)[i])
  
  # datos variable
  x <- diff(df_ts[, i], differences = 1)
  
  # ADF Test 
  adf_type = c("drift", "trend")
  
  for (model_type in adf_type){
    
    adf <- ur.df(x, type = model_type, selectlags = "BIC")
    adf_stat <- adf@teststat[1]
    if (model_type == 'drift') {
      rechazo_1 <- adf_stat < adf@cval[1]
      rechazo_5 <- adf_stat < adf@cval[3]
      rechazo_10 <- adf_stat < adf@cval[5]
      cval1 = adf@cval[1]
      cval5 = adf@cval[3]
      cval10 = adf@cval[5]
    }
    else {
      rechazo_1 <- adf_stat < adf@cval[1]
      rechazo_5 <- adf_stat < adf@cval[4]
      rechazo_10 <- adf_stat < adf@cval[7]
      cval1 = adf@cval[1]
      cval5 = adf@cval[4]
      cval10 = adf@cval[7]
    }
    if (rechazo_1) {
      resultado = sprintf("%f***", adf_stat)
    } 
    else if (rechazo_5) {
      resultado = sprintf("%f**", adf_stat)
    } 
    else if (rechazo_10) {
      resultado = sprintf("%f*", adf_stat)
    } 
    else {
      resultado = sprintf("%f", adf_stat)
    }
    
    #agregamos resultados
    resultados_adf_I1 <- rbind(resultados_adf_I1, data.frame(
      Variable = colnames(df_ts)[i],
      Tipo_Modelo = model_type,
      Estadístico_ADF = resultado,
      cval1 = cval1,
      cval5 = cval5,
      cval10 = cval10,
      stringsAsFactors = FALSE
    ))
    
  }
  
  pp_type = c("constant", "trend")
  # UR PP Test (constant y short)
  
  for (model_type in pp_type){
    pp <- ur.pp(x, type = "Z-tau", model = model_type, lags = 'short')
    stat_pp <- pp@teststat[1]
    
    crit_pp_1 <- pp@cval[1]
    crit_pp_5 <- pp@cval[2]
    crit_pp_10 <- pp@cval[3]
    
    
    rechazo_1 <- stat_pp < crit_pp_1
    rechazo_5 <- stat_pp < crit_pp_5
    rechazo_10 <- stat_pp < crit_pp_10
    if (rechazo_1) {
      resultado_pp = sprintf("%f***", stat_pp)
    } 
    else if (rechazo_5) {
      resultado_pp = sprintf("%f**", stat_pp)
    } 
    else if (rechazo_10) {
      resultado_pp = sprintf("%f*", stat_pp)
    } 
    else {
      resultado_pp = sprintf("%f", stat_pp)
    }
    
    #Agregar resultados a pp
    #agregamos resultados
    resultados_pp_I1 <- rbind(resultados_pp_I1, data.frame(
      Variable = colnames(df_ts)[i],
      Tipo_Modelo = model_type,
      Estadístico_PP = resultado_pp,
      cval1 = pp@cval[1],
      cval5 = pp@cval[2],
      cval10 = pp@cval[3],
      stringsAsFactors = FALSE
    ))
    
    
  }
}


##ZA
#Zivot and Andrews unit root test. Para encontrar breakpoints en las series.
resultados_ZA <- data.frame(
  Variable = character(),
  Tipo_Modelo = character(),
  Estadístico_ZA = character(),
  breakpoin = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(df_ts)) {
  cat("\nVariable:", colnames(df_ts)[i])
  
  # Seleccion de serie
  x <- df_ts[, i]
  
  #tipos de modelos
  za_models = c("intercept", "trend", "both")
  
  # ZA Test 
  for (model_type in za_models){
    za_test <- ur.za(x, model = model_type)
    stat_za <- za_test@teststat
    crit_za_1 <- za_test@cval[1]
    crit_za_5 <- za_test@cval[2]
    crit_za_10 <- za_test@cval[3]
    
    rechazo_1 <- stat_za < crit_za_1
    rechazo_5 <- stat_za < crit_za_5
    rechazo_10 <- stat_za < crit_za_10
    if (rechazo_1) {
      resultado = sprintf("%f***", stat_za)
    } 
    else if (rechazo_5) {
      resultado = sprintf("%f**", stat_za)
    } 
    else if (rechazo_10) {
      resultado = sprintf("%f*", stat_za)
    } 
    else {
      resultado = sprintf("%f", stat_za)
    }
    
    # Encontrar Break point
    
    y_b = za_test@bpoint + 1995 - 1
    
    #Agregar resultados
    #agregamos resultados
    resultados_ZA <- rbind(resultados_ZA, data.frame(
      Variable = colnames(df_ts)[i],
      Tipo_Modelo = model_type,
      Estadístico_ZA = resultado,
      breakpoint = y_b,
      stringsAsFactors = FALSE
    ))
  }
  
  
}


###PROBANDO 

orden_I <- order_integration(df_ts, max_order = 2, method = "boot_ur", level = 0.05, plot_orders = TRUE)




# ====================================
# 4. Pruebas de COINTEGRACION CON ARDL
# ====================================
vars <- c("ECI_SITC","FDI_GDP", "Trade_GDP","GDP_capita","GFCF_GDP",
          "School_enrollment","GovQualityIndex","ICT_Index")
dep      <- "ECI_SITC"           # dependiente fija
fixed    <- "FDI_GDP"            # siempre entra como regresor
cand     <- setdiff(vars, c(dep, fixed))   # candidatos libres
id_counter <- 1
combo_list <- list()
for (c in 1:6) {
  combs <- combn(cand, c, simplify = FALSE)
  for (x_vec in combs) {
    # Guardamos como string ordenada alfabéticamente para asegurar unicidad
    key <- paste(sort(x_vec), collapse = ", ")
    combo_list[[key]] <- id_counter
    id_counter <- id_counter + 1
  }
}

coef_ast_s <- function(coef, p_val) {
   if (p_val < 0.01) {
    return(sprintf("%s***", coef))
  } else if (p_val < 0.05) {
    return(sprintf("%s**", coef))
  } else if (p_val < 0.1) {
    return(sprintf("%s*", coef))
  } else {
    return(sprintf("%s", coef))
  }
}
coef_ast_f <- function(coef, p_val) {
  if (p_val < 0.01) {
    return(sprintf("%f***", coef))
  } else if (p_val < 0.05) {
    return(sprintf("%s**", coef))
  } else if (p_val < 0.1) {
    return(sprintf("%f*", coef))
  } else {
    return(sprintf("%f", coef))
  }
}

safe_extract <- function(expr) {
  tryCatch(expr, error = function(e) "-")
}

# ====================================
# 4.1 Modelo ECI_SITC ~ FDI_GDP
# ====================================


#DATAFRAME
resultados_Cointegracion_ES_F <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  max_lag = character(),
  Estadístico = character(),
  cval1_I0 = numeric(),
  cval5_I0 = numeric(),
  cval10_I0 = numeric(),
  cval1_I1 = numeric(),
  cval5_I1 = numeric(),
  cval10_I1 = numeric(),
  conclusion = character(),
  stringsAsFactors = FALSE
)


result_Longcoef_ES_F <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  Constant = character(),
  FDI = character(),
  TRADE = character(),
  GQI = character(),
  ICT = character(),
  GDPPC = character(),
  SE = character(),
  GFCF = character(),
  Adj_R = character(),
  F_stat = character(),
  AR1 = character(),
  AR2 = character(),
  BP = character(),
  Arch = character(),
  Ramsey_Reset= character(),
  SW_Normality = character(),
  stringsAsFactors = FALSE
)

result_Causality_ES_F <- data.frame(
  id = numeric(),
  Specification = character(),
  Short_Run_FDI = character(),
  Long_Run_ECT = character(),
  stringsAsFactors = FALSE
)
x <- df_ts[, c("ECI_SITC","FDI_GDP", "Trade_GDP","GDP_capita",
               "GFCF_GDP","School_enrollment","GovQualityIndex","ICT_Index")]
vars <- c("ECI_SITC","FDI_GDP", "Trade_GDP","GDP_capita","GFCF_GDP",
          "School_enrollment","GovQualityIndex","ICT_Index")


dep      <- "ECI_SITC"           # dependiente fija
fixed    <- "FDI_GDP"            # siempre entra como regresor
cand     <- setdiff(vars, c(dep, fixed))   # candidatos libres
coef_c <- c("(Intercept)", "ECI_SITC.1", "FDI_GDP.1", "Trade_GDP.1","GDP_capita.1", "GFCF_GDP.1", "School_enrollment.1", "GovQualityIndex.1", "ICT_Index.1")

for (c in 1:6) {
  combs = combn(cand,c)
  for (i in 1:ncol(combs)) {
    x_vec <- combs[,i]
    fml <- as.formula(str_glue("{dep} ~ {fixed} + {paste(sort(x_vec), collapse = ' + ')}"))
    mainvar <- str_glue("{dep} <- {fixed}")
    controls <- str_glue("{paste(sort(x_vec), collapse = ', ')}")
    id_unico <- combo_list[[controls]]
    specif <- str_glue("Model {id_unico}.{dep}/({fixed}, {controls})")
    #Entro al modelo
    if (c == 5 | c == 6) {
      p = 1
      q = 2
    }
    else {
      p = 3
      q = 3
    }
    modelo <- tryCatch({
      ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
    }, error = function(e) {
      message("Error en ardlBound: ", e$message)
      return(NULL)
    })
    if (is.null(modelo)) {
      if (q == 3) {
        p = p - 1
        q = q - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
        }
      else { 
        p = p - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
    }
    
    TTest = FALSE
    while (TTest == FALSE) {
      TTest = TRUE
      
      if(is.nan(bgtest(modelo$model$modelFull$model, order = 1, type = "F")$p.value)){
        TTest = FALSE
      }
      else if(is.nan(bgtest(modelo$model$modelFull$model, order = 2, type = "F")$p.value)){
        TTest = FALSE
      }
      if(is.nan(bptest(modelo$model$modelFull$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(ArchTest(modelo$ARDL.model$residuals, lags = 1)$p.value)){
        TTest = FALSE
      }
      if(is.nan(resettest(modelo$ECM$EC.model$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(modelo$sp$p.value)){
        TTest = FALSE
      }
      
      
      if (TTest == FALSE) {
        
        print("Error en pvalores de tests")
        if (q == 3) {
          p = p - 1
          q = q - 1
        }
        else if (p == 1) {
          q = 1
        }
        
        else { 
          p = p - 1
          }
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
        
      else {
        TTest = TRUE
        }
      }
    
    #Sacamos valores criticos
    out <- capture.output(pssbounds(obs = 28, case = 3, fstat = modelo$F.stat, k = modelo$k))
    linea <- out[12]
    valores_10 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_10_I0 <- valores_10[2]
    cval_10_I1 <- valores_10[3]
    linea <- out[13]
    valores_5 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_5_I0 <- valores_5[2]
    cval_5_I1 <- valores_5[3]
    linea <- out[14]
    valores_1 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_1_I0 <- valores_1[2]
    cval_1_I1 <- valores_1[3]
    
    
    #Conclusion modelo
    if (modelo$F.stat > cval_1_I1) {
      conclusion <- "Cointeg 1%"
      F_stat <- sprintf("%f***", modelo$F.stat) 
    }
    else if (modelo$F.stat > cval_5_I1) {
      conclusion <- "Cointeg 5%"
      F_stat <- sprintf("%f**", modelo$F.stat)
    }
    else if (modelo$F.stat > cval_10_I1) {
      conclusion <- "Cointeg 10%"
      F_stat <- sprintf("%f*", modelo$F.stat)
    }
    else if (modelo$F.stat < cval_10_I0) {
      conclusion <- "No cointeg"
      F_stat <- sprintf("%f", modelo$F.stat)
    }
    else  {
      conclusion <- "Indefinido"
      F_stat <- sprintf("%f", modelo$F.stat)
      
    }
    
    #max lags
    
    lags_vector <- as.numeric(modelo$p)
    lags_str <- sprintf("(%s)", paste(lags_vector, collapse = ","))

    #Resultado
    
    #agregamos resultados Cointegracion
    resultados_Cointegracion_ES_F <- rbind(resultados_Cointegracion_ES_F, data.frame(
      Main = mainvar,
      Control = controls,
      Specification = specif,
      id = id_unico,
      max_lag = lags_str,
      Estadístico = F_stat,
      cval1_I0 = cval_1_I0,
      cval1_I1 = cval_1_I1,
      cval5_I0 = cval_5_I0,
      cval5_I1 = cval_5_I1,
      cval10_I0 = cval_10_I0,
      cval10_I1 = cval_10_I1,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ))
    
    if (conclusion == "Cointeg 10%" | conclusion == "Cointeg 5%" | conclusion == "Cointeg 1%") {
      test_out <- capture.output(summary(modelo$ARDL.model))
      coefs <- test_out[14:25]
      #### SOLUCIONAR HAY VECES QUE COEFS DESDE TEST_OUT NO POSEE LOS DATOS CORRECTOS O NECESARIOS> REC RECORRER TODAS LAS LINEAS Y CHECKEAR
      for (linea in coefs) {
        split <- strsplit(linea, split ="\\s+")
        if (split[[1]][1] %in% coef_c) {
          p_val <- as.numeric(split[[1]][5])
          coef <- split[[1]][2]
          coefficiente <- coef_ast_s(coef, p_val)
          sd <- split[[1]][3]
          if (split[[1]][1] == "(Intercept)") {
            inter_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ECI_SITC.1") {
            ECI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "FDI_GDP.1") {
            FDI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "Trade_GDP.1") {
            Trade_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GovQualityIndex.1") {
            Gov_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ICT_Index.1") {
            ICT_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GDP_capita.1") {
            GDP_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GFCF_GDP.1") {
            GFCF_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "School_enrollment.1") {
            SE_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else {
            print("Error en encontrar coeficiente")
          }}
        
        }
        
        
        
      #F-stat
      F_line <- test_out[length(test_out)-1]
      vals <- as.numeric(unlist(regmatches(F_line, gregexpr("[0-9.]+", F_line))))
      p_val <- vals[4]
      if (p_val < 0.01) {
        F_stat <- sprintf("%f***", vals[1])
      }
      else if (p_val < 0.05) {
        F_stat <- sprintf("%f**", vals[1])
      }
      else if (p_val < 0.1) {
        F_stat <- sprintf("%f*", vals[1])
      }
      else  {
        F_stat <- sprintf("%f", vals[1])
      }
      # R2 Adj
      Rsqr_line <- test_out[length(test_out)-2]
      vals <- as.numeric(unlist(regmatches(Rsqr_line, gregexpr("[0-9.]+", Rsqr_line))))
      r2_adj <- vals[2]
      
      
      #Flag de modelo con algun test no aprobado
      Flag = TRUE
      #AR 1 y 2
      LM_test <- bgtest(modelo$model$modelFull$model, order = 1, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_1 <- coef_ast_f(lmstat, p_val)
      
      AR1 <- str_glue("{AR_1} \n ({p_val})")
      
      LM_test <- bgtest(modelo$model$modelFull$model, order = 2, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_2 <- coef_ast_f(lmstat, p_val)
      
      AR2 <- str_glue("{AR_2} \n ({p_val})")
      
      
      #Breush-Pagan Test for homoskedasticity of residuals
      BP_test <- bptest(modelo$model$modelFull$model)
      p_val <- BP_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      bpstat <- BP_test$statistic
      BP <- coef_ast_f(bpstat, p_val)
      
      BP_st <- str_glue("{BP} \n ({p_val})")
      
      #ARCH test
      Arch_Test <- ArchTest(modelo$ARDL.model$residuals, lags = 1)
      archstat <- Arch_Test$statistic
      p_val <- Arch_Test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      ARCH <- coef_ast_f(archstat, p_val)
      
      ARCH_st <- str_glue("{ARCH} \n ({p_val})")
      
      # Ramsey RESET
      RR_test <- resettest(modelo$ECM$EC.model$model)
      rrstat <- RR_test$statistic
      p_val <- RR_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      RR <- coef_ast_f(rrstat, p_val)
      
      RR_st <- str_glue("{RR} \n ({p_val})")
      
      
      #Shapiro-Wilk test of normality of residuals
      swstat <- modelo$sp$statistic
      p_val <- modelo$sp$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      SW <- coef_ast_f(swstat, p_val)
      
      SW_st <- str_glue("{SW} \n ({p_val})")
      
      
      #GUARDAMOS PLOT CUSUMS
      plot_object <- recordPlot()
      titulo <- paste("\n CUSUMS modelo ", id_unico, " FDI > SITC")
      nombre <- paste("FDI_SITC_",id_unico, ".png")
      png(filename = nombre, width = 1200, height = 800, res = 150)
      replayPlot(plot_object)
      title(titulo)
      dev.off()
      
      ## AQUI DEBERIA GENERAR TABLA
      
      
      result_Longcoef_ES_F <- rbind(result_Longcoef_ES_F, data.frame(
        Main = mainvar,
        Control = controls,
        Specification = specif,
        id = id_unico,
        Constant = inter_coef,
        FDI = safe_extract(FDI_coef),
        TRADE = safe_extract(Trade_coef),
        GQI = safe_extract(Gov_coef),
        ICT = safe_extract(ICT_coef),
        GDPPC = safe_extract(GDP_coef),
        SE = safe_extract(SE_coef),
        GFCF = safe_extract(GFCF_coef),
        Adj_R = r2_adj,
        F_stat = F_stat,
        AR1 = AR1,
        AR2 = AR2,
        BP = BP_st,
        Arch = ARCH_st,
        Ramsey_Reset= RR_st,
        SW_Normality = SW_st,
        stringsAsFactors = FALSE
      ))
        
        
        #AHORA HACER OTRA TABLA
      vars_to_remove <- c("inter_coef", "ECI_coef","FDI_coef", 'Trade_coef',"Gov_coef", "ICT_coef", 
                          "GDP_coef", "SE_coef", "GFCF_coef", "r2_adj", 
                          "F_stat", "AR1", "AR2", "BP_st", "ARCH_st", 
                          "RR_st", "SW_st")
      rm(list = vars_to_remove[vars_to_remove %in% ls()])
      
      
      if (Flag == TRUE) {
        fdis <- c()
        valor <- 0
        for (name in names(coef(modelo$ECM$EC.model))) {
          
          if (grepl("dFDI_GDP",name)) {
            fdis <- c(fdis, name)
            valor <- valor + coef(modelo$ECM$EC.model)[name]
          }
          
        }
        #Short-Run Causality por medio de Joint Significance de dfdis
        src <- linearHypothesis(modelo$ECM$EC.model, paste0(fdis, " = 0"))
        #p valor
        p_value <- src$`Pr(>F)`[2]
        short_run <- coef_ast_f(valor, p_value)
        
        long_run <- coef_ast_f(coef(modelo$ECM$EC.model)["ec.1"],summary(modelo$ECM$EC.model)$coefficients["ec.1", "Pr(>|t|)"])
        result_Causality_ES_F <- rbind(result_Causality_ES_F, data.frame(
          Specification = specif,
          id = id_unico,
          Short_Run_FDI = short_run,
          Long_Run_ECT = long_run,
          stringsAsFactors = FALSE
        ))
      }
      }
    }
}
    


# ====================================
# 4.2 Modelo FDI_GDP ~ ECI_SITC
# ====================================

#DATAFRAME
resultados_Cointegracion_F_ES <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  max_lag = character(),
  Estadístico = character(),
  cval1_I0 = numeric(),
  cval5_I0 = numeric(),
  cval10_I0 = numeric(),
  cval1_I1 = numeric(),
  cval5_I1 = numeric(),
  cval10_I1 = numeric(),
  conclusion = character(),
  stringsAsFactors = FALSE
)


result_Longcoef_F_ES <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  Constant = character(),
  ECI = character(),
  TRADE = character(),
  GQI = character(),
  ICT = character(),
  GDPPC = character(),
  SE = character(),
  GFCF = character(),
  Adj_R = character(),
  F_stat = character(),
  AR1 = character(),
  AR2 = character(),
  BP = character(),
  Arch = character(),
  Ramsey_Reset= character(),
  SW_Normality = character(),
  stringsAsFactors = FALSE
)

result_Causality_F_ES <- data.frame(
  id = numeric(),
  Specification = character(),
  Short_Run_ECI = character(),
  Long_Run_ECT = character(),
  stringsAsFactors = FALSE
)
x <- df_ts[, c("ECI_SITC","FDI_GDP", "Trade_GDP","GDP_capita",
               "GFCF_GDP","School_enrollment","GovQualityIndex","ICT_Index")]
vars <- c("ECI_SITC","FDI_GDP", "Trade_GDP","GDP_capita","GFCF_GDP",
          "School_enrollment","GovQualityIndex","ICT_Index")



fixed      <- "ECI_SITC"           # dependiente fija
dep    <- "FDI_GDP"            # siempre entra como regresor
cand     <- setdiff(vars, c(dep, fixed))   # candidatos libres
coef_c <- c("(Intercept)", "ECI_SITC.1", "FDI_GDP.1", "Trade_GDP.1","GDP_capita.1", "GFCF_GDP.1", "School_enrollment.1", "GovQualityIndex.1", "ICT_Index.1")

for (c in 1:6) {
  combs = combn(cand,c)
  for (i in 1:ncol(combs)) {
    x_vec <- combs[,i]
    fml <- as.formula(str_glue("{dep} ~ {fixed} + {paste(sort(x_vec), collapse = ' + ')}"))
    mainvar <- str_glue("{dep} <- {fixed}")
    controls <- str_glue("{paste(sort(x_vec), collapse = ', ')}")
    id_unico <- combo_list[[controls]]
    specif <- str_glue("Model {id_unico}.{dep}/({fixed}, {controls})")
    #Entro al modelo
    if (c == 5 | c == 6) {
      p = 1
      q = 2
    }
    else {
      p = 3
      q = 3
    }
    modelo <- tryCatch({
      ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
    }, error = function(e) {
      message("Error en ardlBound: ", e$message)
      return(NULL)
    })
    if (is.null(modelo)) {
      if (q == 3) {
        p = p - 1
        q = q - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      else { 
        p = p - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
    }
    
    TTest = FALSE
    while (TTest == FALSE) {
      TTest = TRUE
      
      if(is.nan(bgtest(modelo$model$modelFull$model, order = 1, type = "F")$p.value)){
        TTest = FALSE
      }
      else if(is.nan(bgtest(modelo$model$modelFull$model, order = 2, type = "F")$p.value)){
        TTest = FALSE
      }
      if(is.nan(bptest(modelo$model$modelFull$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(ArchTest(modelo$ARDL.model$residuals, lags = 1)$p.value)){
        TTest = FALSE
      }
      if(is.nan(resettest(modelo$ECM$EC.model$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(modelo$sp$p.value)){
        TTest = FALSE
      }
      
      
      if (TTest == FALSE) {
        
        print("Error en pvalores de tests")
        if (q == 3) {
          p = p - 1
          q = q - 1
        }
        else if (p == 1) {
          q = 1
        }
        else { 
          p = p - 1
        }
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      
      else {
        TTest = TRUE
      }
    }
    
    #Sacamos valores criticos
    out <- capture.output(pssbounds(obs = 28, case = 3, fstat = modelo$F.stat, k = modelo$k))
    linea <- out[12]
    valores_10 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_10_I0 <- valores_10[2]
    cval_10_I1 <- valores_10[3]
    linea <- out[13]
    valores_5 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_5_I0 <- valores_5[2]
    cval_5_I1 <- valores_5[3]
    linea <- out[14]
    valores_1 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_1_I0 <- valores_1[2]
    cval_1_I1 <- valores_1[3]
    
    
    #Conclusion modelo
    if (modelo$F.stat > cval_1_I1) {
      conclusion <- "Cointeg 1%"
      F_stat <- sprintf("%f***", modelo$F.stat) 
    }
    else if (modelo$F.stat > cval_5_I1) {
      conclusion <- "Cointeg 5%"
      F_stat <- sprintf("%f**", modelo$F.stat)
    }
    else if (modelo$F.stat > cval_10_I1) {
      conclusion <- "Cointeg 10%"
      F_stat <- sprintf("%f*", modelo$F.stat)
    }
    else if (modelo$F.stat < cval_10_I0) {
      conclusion <- "No cointeg"
      F_stat <- sprintf("%f", modelo$F.stat)
    }
    else  {
      conclusion <- "Indefinido"
      F_stat <- sprintf("%f", modelo$F.stat)
      
    }
    
    #max lags
    
    lags_vector <- as.numeric(modelo$p)
    lags_str <- sprintf("(%s)", paste(lags_vector, collapse = ","))
    
    #Resultado
    
    #agregamos resultados Cointegracion
    resultados_Cointegracion_F_ES <- rbind(resultados_Cointegracion_F_ES, data.frame(
      Main = mainvar,
      Control = controls,
      Specification = specif,
      id = id_unico,
      max_lag = lags_str,
      Estadístico = F_stat,
      cval1_I0 = cval_1_I0,
      cval1_I1 = cval_1_I1,
      cval5_I0 = cval_5_I0,
      cval5_I1 = cval_5_I1,
      cval10_I0 = cval_10_I0,
      cval10_I1 = cval_10_I1,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ))
    
    if (conclusion == "Cointeg 10%" | conclusion == "Cointeg 5%" | conclusion == "Cointeg 1%") {
      test_out <- capture.output(summary(modelo$ARDL.model))
      coefs <- test_out[14:25]
      #### SOLUCIONAR HAY VECES QUE COEFS DESDE TEST_OUT NO POSEE LOS DATOS CORRECTOS O NECESARIOS> REC RECORRER TODAS LAS LINEAS Y CHECKEAR
      for (linea in coefs) {
        split <- strsplit(linea, split ="\\s+")
        if (split[[1]][1] %in% coef_c) {
          p_val <- as.numeric(split[[1]][5])
          coef <- split[[1]][2]
          coefficiente <- coef_ast_s(coef, p_val)
          sd <- split[[1]][3]
          if (split[[1]][1] == "(Intercept)") {
            inter_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ECI_SITC.1") {
            ECI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "FDI_GDP.1") {
            FDI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "Trade_GDP.1") {
            Trade_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GovQualityIndex.1") {
            Gov_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ICT_Index.1") {
            ICT_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GDP_capita.1") {
            GDP_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GFCF_GDP.1") {
            GFCF_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "School_enrollment.1") {
            SE_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else {
            print("Error en encontrar coeficiente")
          }}
        
      }
      
      
      
      #F-stat
      F_line <- test_out[length(test_out)-1]
      vals <- as.numeric(unlist(regmatches(F_line, gregexpr("[0-9.]+", F_line))))
      p_val <- vals[4]
      if (p_val < 0.01) {
        F_stat <- sprintf("%f***", vals[1])
      }
      else if (p_val < 0.05) {
        F_stat <- sprintf("%f**", vals[1])
      }
      else if (p_val < 0.1) {
        F_stat <- sprintf("%f*", vals[1])
      }
      else  {
        F_stat <- sprintf("%f", vals[1])
      }
      # R2 Adj
      Rsqr_line <- test_out[length(test_out)-2]
      vals <- as.numeric(unlist(regmatches(Rsqr_line, gregexpr("[0-9.]+", Rsqr_line))))
      r2_adj <- vals[2]
      
      
      #Flag de modelo con algun test no aprobado
      Flag = TRUE
      #AR 1 y 2
      LM_test <- bgtest(modelo$model$modelFull$model, order = 1, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_1 <- coef_ast_f(lmstat, p_val)
      
      AR1 <- str_glue("{AR_1} \n ({p_val})")
      
      LM_test <- bgtest(modelo$model$modelFull$model, order = 2, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_2 <- coef_ast_f(lmstat, p_val)
      
      AR2 <- str_glue("{AR_2} \n ({p_val})")
      
      
      #Breush-Pagan Test for homoskedasticity of residuals
      BP_test <- bptest(modelo$model$modelFull$model)
      p_val <- BP_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      bpstat <- BP_test$statistic
      BP <- coef_ast_f(bpstat, p_val)
      
      BP_st <- str_glue("{BP} \n ({p_val})")
      
      #ARCH test
      Arch_Test <- ArchTest(modelo$ARDL.model$residuals, lags = 1)
      archstat <- Arch_Test$statistic
      p_val <- Arch_Test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      ARCH <- coef_ast_f(archstat, p_val)
      
      ARCH_st <- str_glue("{ARCH} \n ({p_val})")
      
      # Ramsey RESET
      RR_test <- resettest(modelo$ECM$EC.model$model)
      rrstat <- RR_test$statistic
      p_val <- RR_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      RR <- coef_ast_f(rrstat, p_val)
      
      RR_st <- str_glue("{RR} \n ({p_val})")
      
      
      #Shapiro-Wilk test of normality of residuals
      swstat <- modelo$sp$statistic
      p_val <- modelo$sp$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      SW <- coef_ast_f(swstat, p_val)
      
      SW_st <- str_glue("{SW} \n ({p_val})")
      
      #GUARDAMOS PLOT CUSUMS
      plot_object <- recordPlot()
      titulo <- paste("\n CUSUMS modelo ", id_unico, "SITC > FDI")
      nombre <- paste("SITC_FDI",id_unico, ".png")
      png(filename = nombre, width = 1200, height = 800, res = 150)
      replayPlot(plot_object)
      title(titulo)
      dev.off()
      
      ## AQUI DEBERIA GENERAR TABLA
      
      
      result_Longcoef_F_ES <- rbind(result_Longcoef_F_ES, data.frame(
        Main = mainvar,
        Control = controls,
        Specification = specif,
        id = id_unico,
        Constant = inter_coef,
        ECI = safe_extract(ECI_coef),
        TRADE = safe_extract(Trade_coef),
        GQI = safe_extract(Gov_coef),
        ICT = safe_extract(ICT_coef),
        GDPPC = safe_extract(GDP_coef),
        SE = safe_extract(SE_coef),
        GFCF = safe_extract(GFCF_coef),
        Adj_R = r2_adj,
        F_stat = F_stat,
        AR1 = AR1,
        AR2 = AR2,
        BP = BP_st,
        Arch = ARCH_st,
        Ramsey_Reset= RR_st,
        SW_Normality = SW_st,
        stringsAsFactors = FALSE
      ))
      
      
      #AHORA HACER OTRA TABLA
      vars_to_remove <- c("inter_coef", "ECI_coef","FDI_coef","Trade_coef", "Gov_coef", "ICT_coef", 
                          "GDP_coef", "SE_coef", "GFCF_coef", "r2_adj", 
                          "F_stat", "AR1", "AR2", "BP_st", "ARCH_st", 
                          "RR_st", "SW_st")
      rm(list = vars_to_remove[vars_to_remove %in% ls()])
      
      
      if (Flag == TRUE) {
        fdis <- c()
        valor <- 0
        for (name in names(coef(modelo$ECM$EC.model))) {
          
          if (grepl("dECI_SITC",name)) {
            fdis <- c(fdis, name)
            valor <- valor + coef(modelo$ECM$EC.model)[name]
          }
          
        }
        #Short-Run Causality por medio de Joint Significance de dfdis
        src <- linearHypothesis(modelo$ECM$EC.model, paste0(fdis, " = 0"))
        #p valor
        p_value <- src$`Pr(>F)`[2]
        short_run <- coef_ast_f(valor, p_value)
        
        long_run <- coef_ast_f(coef(modelo$ECM$EC.model)["ec.1"],summary(modelo$ECM$EC.model)$coefficients["ec.1", "Pr(>|t|)"])
        result_Causality_F_ES <- rbind(result_Causality_F_ES, data.frame(
          Specification = specif,
          id = id_unico,
          Short_Run_ECI = short_run,
          Long_Run_ECT = long_run,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}




# ====================================
# 4.3 Modelo ECI_HS ~ FDI_GDP
# ====================================
#DATAFRAME
resultados_Cointegracion_EH_F <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  max_lag = character(),
  Estadístico = character(),
  cval1_I0 = numeric(),
  cval5_I0 = numeric(),
  cval10_I0 = numeric(),
  cval1_I1 = numeric(),
  cval5_I1 = numeric(),
  cval10_I1 = numeric(),
  conclusion = character(),
  stringsAsFactors = FALSE
)


result_Longcoef_EH_F <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  Constant = character(),
  FDI = character(),
  TRADE = character(),
  GQI = character(),
  ICT = character(),
  GDPPC = character(),
  SE = character(),
  GFCF = character(),
  Adj_R = character(),
  F_stat = character(),
  AR1 = character(),
  AR2 = character(),
  BP = character(),
  Arch = character(),
  Ramsey_Reset= character(),
  SW_Normality = character(),
  stringsAsFactors = FALSE
)

result_Causality_EH_F <- data.frame(
  id = numeric(),
  Specification = character(),
  Short_Run_FDI = character(),
  Long_Run_ECT = character(),
  stringsAsFactors = FALSE
)
x <- df_ts[, c("ECI_HS92","FDI_GDP", "Trade_GDP","GDP_capita",
               "GFCF_GDP","School_enrollment","GovQualityIndex","ICT_Index")]
vars <- c("ECI_HS92","FDI_GDP","Trade_GDP","GDP_capita","GFCF_GDP",
          "School_enrollment","GovQualityIndex","ICT_Index")


dep      <- "ECI_HS92"           # dependiente fija
fixed    <- "FDI_GDP"            # siempre entra como regresor
cand     <- setdiff(vars, c(dep, fixed))   # candidatos libres
coef_c <- c("(Intercept)", "ECI_HS92.1", "FDI_GDP.1", "Trade_GDP.1","GDP_capita.1", "GFCF_GDP.1", "School_enrollment.1", "GovQualityIndex.1", "ICT_Index.1")

for (c in 1:6) {
  combs = combn(cand,c)
  for (i in 1:ncol(combs)) {
    x_vec <- combs[,i]
    fml <- as.formula(str_glue("{dep} ~ {fixed} + {paste(sort(x_vec), collapse = ' + ')}"))
    mainvar <- str_glue("{dep} <- {fixed}")
    controls <- str_glue("{paste(sort(x_vec), collapse = ', ')}")
    id_unico <- combo_list[[controls]]
    specif <- str_glue("Model {id_unico}.{dep}/({fixed}, {controls})")
    #Entro al modelo
    if (c == 5| c == 6) {
      p = 1
      q = 2
    }
    else {
      p = 3
      q = 3
    }
    modelo <- tryCatch({
      ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
    }, error = function(e) {
      message("Error en ardlBound: ", e$message)
      return(NULL)
    })
    if (is.null(modelo)) {
      if (q == 3) {
        p = p - 1
        q = q - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      else { 
        p = p - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
    }
    
    TTest = FALSE
    while (TTest == FALSE) {
      TTest = TRUE
      
      if(is.nan(bgtest(modelo$model$modelFull$model, order = 1, type = "F")$p.value)){
        TTest = FALSE
      }
      else if(is.nan(bgtest(modelo$model$modelFull$model, order = 2, type = "F")$p.value)){
        TTest = FALSE
      }
      if(is.nan(bptest(modelo$model$modelFull$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(ArchTest(modelo$ARDL.model$residuals, lags = 1)$p.value)){
        TTest = FALSE
      }
      if(is.nan(resettest(modelo$ECM$EC.model$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(modelo$sp$p.value)){
        TTest = FALSE
      }
      
      
      if (TTest == FALSE) {
        
        print("Error en pvalores de tests")
        if (q == 3) {
          p = p - 1
          q = q - 1}
        else if (p == 1){
          q = 1
        }
        else { 
          p = p - 1
        }
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      
      else {
        TTest = TRUE
      }
    }
    
    #Sacamos valores criticos
    out <- capture.output(pssbounds(obs = 28, case = 3, fstat = modelo$F.stat, k = modelo$k))
    linea <- out[12]
    valores_10 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_10_I0 <- valores_10[2]
    cval_10_I1 <- valores_10[3]
    linea <- out[13]
    valores_5 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_5_I0 <- valores_5[2]
    cval_5_I1 <- valores_5[3]
    linea <- out[14]
    valores_1 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_1_I0 <- valores_1[2]
    cval_1_I1 <- valores_1[3]
    
    
    #Conclusion modelo
    if (modelo$F.stat > cval_1_I1) {
      conclusion <- "Cointeg 1%"
      F_stat <- sprintf("%f***", modelo$F.stat) 
    }
    else if (modelo$F.stat > cval_5_I1) {
      conclusion <- "Cointeg 5%"
      F_stat <- sprintf("%f**", modelo$F.stat)
    }
    else if (modelo$F.stat > cval_10_I1) {
      conclusion <- "Cointeg 10%"
      F_stat <- sprintf("%f*", modelo$F.stat)
    }
    else if (modelo$F.stat < cval_10_I0) {
      conclusion <- "No cointeg"
      F_stat <- sprintf("%f", modelo$F.stat)
    }
    else  {
      conclusion <- "Indefinido"
      F_stat <- sprintf("%f", modelo$F.stat)
      
    }
    
    #max lags
    
    lags_vector <- as.numeric(modelo$p)
    lags_str <- sprintf("(%s)", paste(lags_vector, collapse = ","))
    
    #Resultado
    
    #agregamos resultados Cointegracion
    resultados_Cointegracion_EH_F <- rbind(resultados_Cointegracion_EH_F, data.frame(
      Main = mainvar,
      Control = controls,
      Specification = specif,
      id = id_unico,
      max_lag = lags_str,
      Estadístico = F_stat,
      cval1_I0 = cval_1_I0,
      cval1_I1 = cval_1_I1,
      cval5_I0 = cval_5_I0,
      cval5_I1 = cval_5_I1,
      cval10_I0 = cval_10_I0,
      cval10_I1 = cval_10_I1,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ))
    
    if (conclusion == "Cointeg 10%" | conclusion == "Cointeg 5%" | conclusion == "Cointeg 1%") {
      test_out <- capture.output(summary(modelo$ARDL.model))
      coefs <- test_out[14:25]
      #### SOLUCIONAR HAY VECES QUE COEFS DESDE TEST_OUT NO POSEE LOS DATOS CORRECTOS O NECESARIOS> REC RECORRER TODAS LAS LINEAS Y CHECKEAR
      for (linea in coefs) {
        split <- strsplit(linea, split ="\\s+")
        if (split[[1]][1] %in% coef_c) {
          p_val <- as.numeric(split[[1]][5])
          coef <- split[[1]][2]
          coefficiente <- coef_ast_s(coef, p_val)
          sd <- split[[1]][3]
          if (split[[1]][1] == "(Intercept)") {
            inter_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ECI_HS92.1") {
            ECI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "FDI_GDP.1") {
            FDI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "Trade_GDP.1") {
            Trade_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GovQualityIndex.1") {
            Gov_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ICT_Index.1") {
            ICT_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GDP_capita.1") {
            GDP_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GFCF_GDP.1") {
            GFCF_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "School_enrollment.1") {
            SE_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else {
            print("Error en encontrar coeficiente")
          }}
        
      }
      
      
      
      #F-stat
      F_line <- test_out[length(test_out)-1]
      vals <- as.numeric(unlist(regmatches(F_line, gregexpr("[0-9.]+", F_line))))
      p_val <- vals[4]
      if (p_val < 0.01) {
        F_stat <- sprintf("%f***", vals[1])
      }
      else if (p_val < 0.05) {
        F_stat <- sprintf("%f**", vals[1])
      }
      else if (p_val < 0.1) {
        F_stat <- sprintf("%f*", vals[1])
      }
      else  {
        F_stat <- sprintf("%f", vals[1])
      }
      # R2 Adj
      Rsqr_line <- test_out[length(test_out)-2]
      vals <- as.numeric(unlist(regmatches(Rsqr_line, gregexpr("[0-9.]+", Rsqr_line))))
      r2_adj <- vals[2]
      
      
      #Flag de modelo con algun test no aprobado
      Flag = TRUE
      #AR 1 y 2
      LM_test <- bgtest(modelo$model$modelFull$model, order = 1, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_1 <- coef_ast_f(lmstat, p_val)
      
      AR1 <- str_glue("{AR_1} \n ({p_val})")
      
      LM_test <- bgtest(modelo$model$modelFull$model, order = 2, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_2 <- coef_ast_f(lmstat, p_val)
      
      AR2 <- str_glue("{AR_2} \n ({p_val})")
      
      
      #Breush-Pagan Test for homoskedasticity of residuals
      BP_test <- bptest(modelo$model$modelFull$model)
      p_val <- BP_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      bpstat <- BP_test$statistic
      BP <- coef_ast_f(bpstat, p_val)
      
      BP_st <- str_glue("{BP} \n ({p_val})")
      
      #ARCH test
      Arch_Test <- ArchTest(modelo$ARDL.model$residuals, lags = 1)
      archstat <- Arch_Test$statistic
      p_val <- Arch_Test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      ARCH <- coef_ast_f(archstat, p_val)
      
      ARCH_st <- str_glue("{ARCH} \n ({p_val})")
      
      # Ramsey RESET
      RR_test <- resettest(modelo$ECM$EC.model$model)
      rrstat <- RR_test$statistic
      p_val <- RR_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      RR <- coef_ast_f(rrstat, p_val)
      
      RR_st <- str_glue("{RR} \n ({p_val})")
      
      
      #Shapiro-Wilk test of normality of residuals
      swstat <- modelo$sp$statistic
      p_val <- modelo$sp$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      SW <- coef_ast_f(swstat, p_val)
      
      SW_st <- str_glue("{SW} \n ({p_val})")
      
      
      #GUARDAMOS PLOT CUSUMS
      plot_object <- recordPlot()
      titulo <- paste("\n CUSUMS modelo ", id_unico, "FDI > HS92")
      nombre <- paste("FDI_HS_",id_unico, ".png")
      png(filename = nombre, width = 1200, height = 800, res = 150)
      replayPlot(plot_object)
      title(titulo)
      dev.off()
      ## AQUI DEBERIA GENERAR TABLA
      
      
      result_Longcoef_EH_F <- rbind(result_Longcoef_EH_F, data.frame(
        Main = mainvar,
        Control = controls,
        Specification = specif,
        id = id_unico,
        Constant = inter_coef,
        FDI = safe_extract(FDI_coef),
        TRADE = safe_extract(Trade_coef),
        GQI = safe_extract(Gov_coef),
        ICT = safe_extract(ICT_coef),
        GDPPC = safe_extract(GDP_coef),
        SE = safe_extract(SE_coef),
        GFCF = safe_extract(GFCF_coef),
        Adj_R = r2_adj,
        F_stat = F_stat,
        AR1 = AR1,
        AR2 = AR2,
        BP = BP_st,
        Arch = ARCH_st,
        Ramsey_Reset= RR_st,
        SW_Normality = SW_st,
        stringsAsFactors = FALSE
      ))
      
      
      #AHORA HACER OTRA TABLA
      vars_to_remove <- c("inter_coef", "ECI_coef","FDI_coef", "Trade_coef","Gov_coef", "ICT_coef", 
                          "GDP_coef", "SE_coef", "GFCF_coef", "r2_adj", 
                          "F_stat", "AR1", "AR2", "BP_st", "ARCH_st", 
                          "RR_st", "SW_st")
      rm(list = vars_to_remove[vars_to_remove %in% ls()])
      
      
      if (Flag == TRUE) {
        fdis <- c()
        valor <- 0
        for (name in names(coef(modelo$ECM$EC.model))) {
          
          if (grepl("dFDI_GDP",name)) {
            fdis <- c(fdis, name)
            valor <- valor + coef(modelo$ECM$EC.model)[name]
          }
          
        }
        #Short-Run Causality por medio de Joint Significance de dfdis
        src <- linearHypothesis(modelo$ECM$EC.model, paste0(fdis, " = 0"))
        #p valor
        p_value <- src$`Pr(>F)`[2]
        short_run <- coef_ast_f(valor, p_value)
        
        long_run <- coef_ast_f(coef(modelo$ECM$EC.model)["ec.1"],summary(modelo$ECM$EC.model)$coefficients["ec.1", "Pr(>|t|)"])
        result_Causality_EH_F <- rbind(result_Causality_EH_F, data.frame(
          Specification = specif,
          id = id_unico,
          Short_Run_FDI = short_run,
          Long_Run_ECT = long_run,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# ====================================
# 4.4 Modelo FDI_GDP ~ ECI_HS
# ====================================
#DATAFRAME
resultados_Cointegracion_F_EH <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  max_lag = character(),
  Estadístico = character(),
  cval1_I0 = numeric(),
  cval5_I0 = numeric(),
  cval10_I0 = numeric(),
  cval1_I1 = numeric(),
  cval5_I1 = numeric(),
  cval10_I1 = numeric(),
  conclusion = character(),
  stringsAsFactors = FALSE
)


result_Longcoef_F_EH <- data.frame(
  Main = character(),
  Control = character(),
  id = numeric(),
  Specification = character(),
  Constant = character(),
  ECI = character(),
  TRADE = character(),
  GQI = character(),
  ICT = character(),
  GDPPC = character(),
  SE = character(),
  GFCF = character(),
  Adj_R = character(),
  F_stat = character(),
  AR1 = character(),
  AR2 = character(),
  BP = character(),
  Arch = character(),
  Ramsey_Reset= character(),
  SW_Normality = character(),
  stringsAsFactors = FALSE
)

result_Causality_F_EH <- data.frame(
  id = numeric(),
  Specification = character(),
  Short_Run_ECI = character(),
  Long_Run_ECT = character(),
  stringsAsFactors = FALSE
)
x <- df_ts[, c("ECI_HS92","FDI_GDP", "Trade_GDP","GDP_capita",
               "GFCF_GDP","School_enrollment","GovQualityIndex","ICT_Index")]
vars <- c("ECI_HS92","FDI_GDP", "Trade_GDP","GDP_capita","GFCF_GDP",
          "School_enrollment","GovQualityIndex","ICT_Index")



fixed      <- "ECI_HS92"           # dependiente fija
dep    <- "FDI_GDP"            # siempre entra como regresor
cand     <- setdiff(vars, c(dep, fixed))   # candidatos libres
coef_c <- c("(Intercept)", "ECI_HS92.1", "FDI_GDP.1", "Trade_GDP.1","GDP_capita.1", "GFCF_GDP.1", "School_enrollment.1", "GovQualityIndex.1", "ICT_Index.1")

for (c in 1:6) {
  combs = combn(cand,c)
  for (i in 1:ncol(combs)) {
    x_vec <- combs[,i]
    fml <- as.formula(str_glue("{dep} ~ {fixed} + {paste(sort(x_vec), collapse = ' + ')}"))
    mainvar <- str_glue("{dep} <- {fixed}")
    controls <- str_glue("{paste(sort(x_vec), collapse = ', ')}")
    id_unico <- combo_list[[controls]]
    specif <- str_glue("Model {id_unico}.{dep}/({fixed}, {controls})")
    #Entro al modelo
    if (c == 5| c == 6) {
      p = 1
      q = 2
    }
    else {
      p = 3
      q = 3
    }
    modelo <- tryCatch({
      ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
    }, error = function(e) {
      message("Error en ardlBound: ", e$message)
      return(NULL)
    })
    if (is.null(modelo)) {
      if (q == 3) {
        p = p - 1
        q = q - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      else { 
        p = p - 1
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
    }
    
    TTest = FALSE
    while (TTest == FALSE) {
      TTest = TRUE
      
      if(is.nan(bgtest(modelo$model$modelFull$model, order = 1, type = "F")$p.value)){
        TTest = FALSE
      }
      else if(is.nan(bgtest(modelo$model$modelFull$model, order = 2, type = "F")$p.value)){
        TTest = FALSE
      }
      if(is.nan(bptest(modelo$model$modelFull$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(ArchTest(modelo$ARDL.model$residuals, lags = 1)$p.value)){
        TTest = FALSE
      }
      if(is.nan(resettest(modelo$ECM$EC.model$model)$p.value)){
        TTest = FALSE
      }
      if(is.nan(modelo$sp$p.value)){
        TTest = FALSE
      }
      
      
      if (TTest == FALSE) {
        
        print("Error en pvalores de tests")
        if (q == 3) {
          p = p - 1
          q = q - 1}
        else if (p == 1) {
          q = 1
        }
        else { 
          p = p - 1
        }
        modelo <- ardlBound(data = x, formula = fml, autoOrder = TRUE, ic = "BIC", max.p = p, max.q = q, case = 3, ECM = TRUE, HAC = TRUE, stability = TRUE)
      }
      
      else {
        TTest = TRUE
      }
    }
    
    #Sacamos valores criticos
    out <- capture.output(pssbounds(obs = 28, case = 3, fstat = modelo$F.stat, k = modelo$k))
    linea <- out[12]
    valores_10 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_10_I0 <- valores_10[2]
    cval_10_I1 <- valores_10[3]
    linea <- out[13]
    valores_5 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_5_I0 <- valores_5[2]
    cval_5_I1 <- valores_5[3]
    linea <- out[14]
    valores_1 <- as.numeric(unlist(regmatches(linea, gregexpr("[0-9.]+", linea))))
    cval_1_I0 <- valores_1[2]
    cval_1_I1 <- valores_1[3]
    
    
    #Conclusion modelo
    if (modelo$F.stat > cval_1_I1) {
      conclusion <- "Cointeg 1%"
      F_stat <- sprintf("%f***", modelo$F.stat) 
    }
    else if (modelo$F.stat > cval_5_I1) {
      conclusion <- "Cointeg 5%"
      F_stat <- sprintf("%f**", modelo$F.stat)
    }
    else if (modelo$F.stat > cval_10_I1) {
      conclusion <- "Cointeg 10%"
      F_stat <- sprintf("%f*", modelo$F.stat)
    }
    else if (modelo$F.stat < cval_10_I0) {
      conclusion <- "No cointeg"
      F_stat <- sprintf("%f", modelo$F.stat)
    }
    else  {
      conclusion <- "Indefinido"
      F_stat <- sprintf("%f", modelo$F.stat)
      
    }
    
    #max lags
    
    lags_vector <- as.numeric(modelo$p)
    lags_str <- sprintf("(%s)", paste(lags_vector, collapse = ","))
    
    #Resultado
    
    #agregamos resultados Cointegracion
    resultados_Cointegracion_F_EH <- rbind(resultados_Cointegracion_F_EH, data.frame(
      Main = mainvar,
      Control = controls,
      Specification = specif,
      id = id_unico,
      max_lag = lags_str,
      Estadístico = F_stat,
      cval1_I0 = cval_1_I0,
      cval1_I1 = cval_1_I1,
      cval5_I0 = cval_5_I0,
      cval5_I1 = cval_5_I1,
      cval10_I0 = cval_10_I0,
      cval10_I1 = cval_10_I1,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ))
    
    if (conclusion == "Cointeg 10%" | conclusion == "Cointeg 5%" | conclusion == "Cointeg 1%") {
      test_out <- capture.output(summary(modelo$ARDL.model))
      coefs <- test_out[14:25]
      #### SOLUCIONAR HAY VECES QUE COEFS DESDE TEST_OUT NO POSEE LOS DATOS CORRECTOS O NECESARIOS> REC RECORRER TODAS LAS LINEAS Y CHECKEAR
      for (linea in coefs) {
        split <- strsplit(linea, split ="\\s+")
        if (split[[1]][1] %in% coef_c) {
          p_val <- as.numeric(split[[1]][5])
          coef <- split[[1]][2]
          coefficiente <- coef_ast_s(coef, p_val)
          sd <- split[[1]][3]
          if (split[[1]][1] == "(Intercept)") {
            inter_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ECI_HS92.1") {
            ECI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "FDI_GDP.1") {
            FDI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "Trade_GDP.1") {
            FDI_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GovQualityIndex.1") {
            Gov_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "ICT_Index.1") {
            ICT_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GDP_capita.1") {
            GDP_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "GFCF_GDP.1") {
            GFCF_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else if (split[[1]][1] == "School_enrollment.1") {
            SE_coef <- str_glue("{coefficiente} \n ({sd})")
          }
          else {
            print("Error en encontrar coeficiente")
          }}
        
      }
      
      
      
      #F-stat
      F_line <- test_out[length(test_out)-1]
      vals <- as.numeric(unlist(regmatches(F_line, gregexpr("[0-9.]+", F_line))))
      p_val <- vals[4]
      if (p_val < 0.01) {
        F_stat <- sprintf("%f***", vals[1])
      }
      else if (p_val < 0.05) {
        F_stat <- sprintf("%f**", vals[1])
      }
      else if (p_val < 0.1) {
        F_stat <- sprintf("%f*", vals[1])
      }
      else  {
        F_stat <- sprintf("%f", vals[1])
      }
      # R2 Adj
      Rsqr_line <- test_out[length(test_out)-2]
      vals <- as.numeric(unlist(regmatches(Rsqr_line, gregexpr("[0-9.]+", Rsqr_line))))
      r2_adj <- vals[2]
      
      
      #Flag de modelo con algun test no aprobado
      Flag = TRUE
      #AR 1 y 2
      LM_test <- bgtest(modelo$model$modelFull$model, order = 1, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_1 <- coef_ast_f(lmstat, p_val)
      
      AR1 <- str_glue("{AR_1} \n ({p_val})")
      
      LM_test <- bgtest(modelo$model$modelFull$model, order = 2, type = "F")
      p_val <- LM_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      lmstat <- LM_test$statistic
      AR_2 <- coef_ast_f(lmstat, p_val)
      
      AR2 <- str_glue("{AR_2} \n ({p_val})")
      
      
      #Breush-Pagan Test for homoskedasticity of residuals
      BP_test <- bptest(modelo$model$modelFull$model)
      p_val <- BP_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      bpstat <- BP_test$statistic
      BP <- coef_ast_f(bpstat, p_val)
      
      BP_st <- str_glue("{BP} \n ({p_val})")
      
      #ARCH test
      Arch_Test <- ArchTest(modelo$ARDL.model$residuals, lags = 1)
      archstat <- Arch_Test$statistic
      p_val <- Arch_Test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      ARCH <- coef_ast_f(archstat, p_val)
      
      ARCH_st <- str_glue("{ARCH} \n ({p_val})")
      
      # Ramsey RESET
      RR_test <- resettest(modelo$ECM$EC.model$model)
      rrstat <- RR_test$statistic
      p_val <- RR_test$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      RR <- coef_ast_f(rrstat, p_val)
      
      RR_st <- str_glue("{RR} \n ({p_val})")
      
      
      #Shapiro-Wilk test of normality of residuals
      swstat <- modelo$sp$statistic
      p_val <- modelo$sp$p.value
      if (p_val<0.1) {
        Flag = FALSE
      }
      SW <- coef_ast_f(swstat, p_val)
      
      SW_st <- str_glue("{SW} \n ({p_val})")
      
      #GUARDAMOS PLOT CUSUMS
      plot_object <- recordPlot()
      titulo <- paste("\n CUSUMS modelo ", id_unico, "HS92 > FDI")
      nombre <- paste("HS_FDI",id_unico, ".png")
      png(filename = nombre, width = 1200, height = 800, res = 150)
      replayPlot(plot_object)
      title(titulo)
      dev.off()
      
      ## AQUI DEBERIA GENERAR TABLA
      
      
      result_Longcoef_F_EH <- rbind(result_Longcoef_F_EH, data.frame(
        Main = mainvar,
        Control = controls,
        Specification = specif,
        id = id_unico,
        Constant = inter_coef,
        ECI = safe_extract(ECI_coef),
        TRADE = safe_extract(Trade_coef),
        GQI = safe_extract(Gov_coef),
        ICT = safe_extract(ICT_coef),
        GDPPC = safe_extract(GDP_coef),
        SE = safe_extract(SE_coef),
        GFCF = safe_extract(GFCF_coef),
        Adj_R = r2_adj,
        F_stat = F_stat,
        AR1 = AR1,
        AR2 = AR2,
        BP = BP_st,
        Arch = ARCH_st,
        Ramsey_Reset= RR_st,
        SW_Normality = SW_st,
        stringsAsFactors = FALSE
      ))
      
      
      #AHORA HACER OTRA TABLA
      vars_to_remove <- c("inter_coef", "ECI_coef","FDI_coef", "Trade_coef", "Gov_coef", "ICT_coef", 
                          "GDP_coef", "SE_coef", "GFCF_coef", "r2_adj", 
                          "F_stat", "AR1", "AR2", "BP_st", "ARCH_st", 
                          "RR_st", "SW_st")
      rm(list = vars_to_remove[vars_to_remove %in% ls()])
      
      
      if (Flag == TRUE) {
        fdis <- c()
        valor <- 0
        for (name in names(coef(modelo$ECM$EC.model))) {
          
          if (grepl("dECI_HS92",name)) {
            fdis <- c(fdis, name)
            valor <- valor + coef(modelo$ECM$EC.model)[name]
          }
          
        }
        #Short-Run Causality por medio de Joint Significance de dfdis
        src <- linearHypothesis(modelo$ECM$EC.model, paste0(fdis, " = 0"))
        #p valor
        p_value <- src$`Pr(>F)`[2]
        short_run <- coef_ast_f(valor, p_value)
        
        long_run <- coef_ast_f(coef(modelo$ECM$EC.model)["ec.1"],summary(modelo$ECM$EC.model)$coefficients["ec.1", "Pr(>|t|)"])
        result_Causality_F_EH <- rbind(result_Causality_F_EH, data.frame(
          Specification = specif,
          id = id_unico,
          Short_Run_ECI = short_run,
          Long_Run_ECT = long_run,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}








# ====================================
# 5 utils
# ====================================
# install.packages("xtable")
library(readr)
write_csv(resultados_Cointegracion_F_ES, "C:/Users/bruno/Desktop/Python Tesis/coint2.csv")


write_csv(result_Longcoef_ES_F, "C:/Users/bruno/Desktop/Python Tesis/Longcoef1.csv")
write_csv(result_Longcoef_F_ES, "C:/Users/bruno/Desktop/Python Tesis/Longcoef2.csv")
write_csv(result_Longcoef_EH_F, "C:/Users/bruno/Desktop/Python Tesis/Longcoef3.csv")
write_csv(result_Longcoef_F_EH, "C:/Users/bruno/Desktop/Python Tesis/Longcoef4.csv")


write_csv(result_Causality_ES_F, "C:/Users/bruno/Desktop/Python Tesis/caus1.csv")
write_csv(result_Causality_F_ES, "C:/Users/bruno/Desktop/Python Tesis/caus2.csv")
write_csv(result_Causality_EH_F, "C:/Users/bruno/Desktop/Python Tesis/caus3.csv")
write_csv(result_Causality_F_EH, "C:/Users/bruno/Desktop/Python Tesis/caus4.csv")
