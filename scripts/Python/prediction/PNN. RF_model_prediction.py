# Cargar bibliotecas de Python
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.svm import SVC
import pandas as pd
import numpy as np
import joblib

# Cargar modelo entrenado
filename = "E:/ME/ALL/Project_NN_proteins/SCRIPTS/PYTHON/MODELS/idr_rf_model.sav"
idr_model = joblib.load(filename)

# Cargar datos de la nueva IDR a predecir
df_IDRs = pd.read_csv("E:/ME/ALL/Project_NN_proteins/DATA/DE_NOVO/GOOSE/DATA/CLEAN/IDRBS_DN_goose_kappa_sparrow_890.csv")

# Obtener los datos numéricos para entrenar el modelo SVM
pred_IDR = df_IDRs.iloc[:,[3, 4, 5]].values

# Definir un "scaler"
scaler = StandardScaler()

# Estandarizar los datos para la predicción
pred_IDR = scaler.fit_transform(pred_IDR)

# Realizar predicción de la IDR
pred_result = idr_model.predict(pred_IDR)

# Convertir resultados en data frame
df_pred = pd.DataFrame(pred_result.tolist(), columns = ["Prediction"])

# Juntar data frames
df_IDRs.insert(len(df_IDRs.columns), "Prediction", df_pred)

# Remover columna no asignada
#df_IDRs = df_IDRs.drop(columns=["Unnamed: 0"])

# Guardar data frame
df_IDRs.to_csv("E:/ME/ALL/Project_NN_proteins/DATA/DE_NOVO/SVM/IDRBS_DN_RF_goose_kappa_890.csv")