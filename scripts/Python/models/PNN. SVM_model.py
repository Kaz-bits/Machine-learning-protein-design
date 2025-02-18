# Importar bibliotecas
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report
from sklearn.metrics import RocCurveDisplay
from sklearn.datasets import load_digits
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import joblib
import random

# Cargar archivo de la biblioteca de las IDRs
df_delta = pd.read_csv("D:/MASTER_FILES/DATA/ML/IDRBS_library_sparrow_141_NN_SVM1.csv")

# Grid search ----
# Obtener los datos numéricos para entrenar el modelo SVM
x = df_delta.iloc[:,[3, 4, 5]].values

# Obtener la variable de clasificación
y = df_delta.iloc[:, 2].values

# Crear el pipeline para escalar y ajustar el modelo
pipeline = make_pipeline(StandardScaler(), SVC(kernel='rbf'))

# Definir el rango de valores para C y gamma
param_grid = {
    'svc__C': [0.1, 1, 10, 100, 1000],
    'svc__gamma': [0.001, 0.01, 0.1, 1, 10]
}

# Configurar la búsqueda en cuadrícula con validación cruzada
grid_search = GridSearchCV(pipeline, param_grid, cv=7, scoring='accuracy', n_jobs=-1)

# Ajustar el modelo y encontrar los mejores parámetros
grid_search.fit(x, y)

# Obtener los mejores parámetros y la mejor precisión
print("Mejores parámetros:", grid_search.best_params_)
print("Mejor precisión:", grid_search.best_score_)




# Training ----
# Configuración de la semilla
# temp = random.randint(1, 10000)

# Select seed
temp = 5828

# Dividir los datos en conjuntos de entrenamiento y prueba
x_train, x_test, y_train, y_test = train_test_split(x, y, shuffle = True, test_size = 0.20, random_state = temp)

# Defninir un "scaler"
scaler = StandardScaler()

# Escalar los datos de entrenamiento
x_train = scaler.fit_transform(x_train)

# Escalar los datos de prueba
x_test = scaler.transform(x_test)

# Crear el modelo SVM
model = SVC(C = 1, gamma = 0.01, kernel='rbf')

# Entrenar el modelo con los datos de entrenamiento
model.fit(x_train, y_train)

# Realizar predicciones en el conjunto de prueba
y_pred = model.predict(x_test)

# Calcular la exactitud del modelo
accuracy = accuracy_score(y_test, y_pred)
print("Exactitud del modelo:", round(accuracy, 4))

# Calcular la precision del modelo 
precision = precision_score(y_test, y_pred, pos_label = "High")
print("Precisión del modelo:", round(precision, 4))

# Determinar el recall (true positives)
recall = recall_score(y_test, y_pred, pos_label = "High")
print("Verdaderos positivos:", round(recall, 4))

# Calcular el valor de f1
score_f1 = f1_score(y_test, y_pred, pos_label = "High")
print("F1 score:", round(score_f1, 4))

# Mostrar el reporte de clasificación del modelo SVM
print(classification_report(y_test, y_pred))

# Calcular el valor del área bajo la curva del gráfico ROC (AUC)
svc = RocCurveDisplay.from_estimator(model, x_test, y_test, pos_label = "High")
plt.show()


# ======================= Validación cruzada ==========================

# Definir el número de folds
num_folds = 7

# Crear el generador de validación cruzada k-fold estratificada
cv_stratified = StratifiedKFold(n_splits = num_folds)

# Inicializar lista para almacenar los resultados de las métricas
accuracy_scores = []
precision_scores = []
recall_scores = []
f1_scores = []

# Realizar la validación cruzada
for train_index, test_index in cv_stratified.split(x_train, y_train):
    x_train_fold, x_test_fold = x_train[train_index], x_train[test_index]
    y_train_fold, y_test_fold = y_train[train_index], y_train[test_index]

    # Escalar los datos de entrenamiento y validación de cada fold
    scaler_fold = StandardScaler()
    x_train_fold = scaler_fold.fit_transform(x_train_fold)
    x_test_fold = scaler_fold.transform(x_test_fold)

    # Train the model in the current fold
    model.fit(x_train_fold, y_train_fold)

    # Realizar predicciones en el conjunto de validación
    y_pred_fold = model.predict(x_test_fold)

    # Calcular la precisión del modelo en el fold actual
    accuracy_scores.append(accuracy_score(y_test_fold, y_pred_fold))
    precision_scores.append(precision_score(y_test_fold, y_pred_fold, pos_label = "High"))
    recall_scores.append(recall_score(y_test_fold, y_pred_fold, pos_label = "High"))
    f1_scores.append(f1_score(y_test_fold, y_pred_fold, pos_label = "High"))

# Calcular la precisión promedio y desviación estándar de las métricas de cada fold
# Calcular promedio y desviación estándar para cada métrica
mean_accuracy = np.mean(accuracy_scores)
std_accuracy = np.std(accuracy_scores)
mean_precision = np.mean(precision_scores)
std_precision = np.std(precision_scores)
mean_recall = np.mean(recall_scores)
std_recall = np.std(recall_scores)
mean_f1 = np.mean(f1_scores)
std_f1 = np.std(f1_scores)

# Imprimir los resultados de validación cruzada
print("Resultados de validación cruzada:")
print("Exactitud promedio:", round(mean_accuracy, 2), "±", round(std_accuracy, 2))
print("Precisión promedio:", round(mean_precision, 2), "±", round(std_precision, 2))
print("Recall promedio:", round(mean_recall, 2), "±", round(std_recall, 2))
print("F1 score promedio:", round(mean_f1, 2), "±", round(std_f1, 2))

# Print seed
print("The seed for this model is:", temp)

# Guardar modelo entrenado
filename = "D:/MASTER_FILES/SCRIPTS/PYTHON/MODELS/idr_svm_model.sav"
joblib.dump(model, filename)