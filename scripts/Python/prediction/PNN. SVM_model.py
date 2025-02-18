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
from sklearn.model_selection import ValidationCurveDisplay
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import joblib

# Cargar archivo de la biblioteca de las IDRs
df_delta = pd.read_csv("E:/ME/ALL/Project_NN_proteins/DATA/DATABASES/IDRBS_library_sparrow_141_NN_SVM.csv")

# Obtener los datos numéricos para entrenar el modelo SVM
x = df_delta.iloc[:,[4, 5, 16]].values

# Obtener la variable de clasificación
y = df_delta.iloc[:, 2].values

# Dividir los datos en conjuntos de entrenamiento y prueba
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.20, random_state = 17)

# Defninir un "scaler"
scaler = StandardScaler()

# Escalar los datos de entrenamiento
x_train = scaler.fit_transform(x_train)

# Escalar los datos de prueba
x_test = scaler.fit_transform(x_test)

# Crear el modelo SVM
model = SVC(C = 1, gamma = 0.01, kernel='rbf')

# Entrenar el modelo con los datos de entrenamiento
model.fit(x_train, y_train)

# Realizar predicciones en el conjunto de prueba
y_pred = model.predict(x_test)

# Calcular la exactitud del modelo
accuracy = accuracy_score(y_test, y_pred)
print("Exactitud del modelo:", round(accuracy, 2))

# Calcular la precision del modelo 
precision = precision_score(y_test, y_pred, pos_label = "Alta")
print("Precisión del modelo:", round(precision, 2))

# Determinar el recall (true positives)
recall = recall_score(y_test, y_pred, pos_label = "Alta")
print("Verdaderos positivos:", round(recall, 2))

# Calcular el valor de f1
score_f1 = f1_score(y_test, y_pred, pos_label = "Alta")
print("F1 score:", round(score_f1, 2))


# ======================= Validación cruzada ==========================

# Crear clasificador SVM
svm_classifier = SVC(C = 1, gamma = 0.01, kernel='rbf')

# Definir el número de folds
num_folds = 5

# Crear el generador de validación cruzada k-fold estratificada
cv_stratified = StratifiedKFold(n_splits = num_folds, shuffle = True, random_state = 17)

# Inicializar lista para almacenar los resultados de precisión
accuracy_scores = []

# Realizar la validación cruzada
for train_index, test_index in cv_stratified.split(x_train, y_train):
    x_train_fold, x_val_fold = x_train[train_index], x_train[test_index]
    y_train_fold, y_val_fold = y_train[train_index], y_train[test_index]

    # Entrenar el modelo SVM en el fold actual
    svm_classifier.fit(x_train_fold, y_train_fold)

    # Realizar predicciones en el conjunto de validación
    y_pred_fold = svm_classifier.predict(x_val_fold)

    # Calcular la precisión del modelo en el fold actual
    accuracy = accuracy_score(y_val_fold, y_pred_fold)
    accuracy_scores.append(accuracy)

# Calcular la precisión promedio y desviación estándar de las fold
mean_accuracy = np.mean(accuracy_scores)
std_accuracy = np.std(accuracy_scores)

# Imprimir los resultados
print("Precisión promedio:", round(mean_accuracy, 3))
print("Desviación estándar de la precisión:", round(std_accuracy, 2))

# Mostrar el reporte de clasificación del modelo SVM
print(classification_report(y_test, y_pred))

# Calcular el valor del área bajo la curva del gráfico ROC (AUC)
svc = RocCurveDisplay.from_estimator(model, x_test, y_test, pos_label = "Alta")
plt.show()

# Guardar modelo entrenado
filename = "E:/ME/ALL/Project_NN_proteins/SCRIPTS/PYTHON/MODELS/idr_svm_model.sav"
joblib.dump(model, filename)