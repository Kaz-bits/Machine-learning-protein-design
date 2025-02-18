# Data Processing
import pandas as pd
import numpy as np

# Modelling
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint

# Tree Visualisation
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz

# Cargar archivo de la biblioteca de las IDRs
df_delta = pd.read_csv("E:/ME/ALL/Project_NN_proteins/DATA/DATABASES/IDRBS_library_sparrow_141_NN_SVM.csv")

# Obtener los datos numéricos para entrenar el modelo SVM
x = df_delta.iloc[:,[4, 5, 16]].values

# Obtener la variable de clasificación
y = df_delta.iloc[:, 2].values

# Dividir los datos en conjuntos de entrenamiento y prueba
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.20, random_state = 1)

# Entrenar el modelo
rf = RandomForestClassifier(n_estimators = 222, max_depth = 17)
rf.fit(x_train, y_train)

# Realizar predicciones
y_pred = rf.predict(x_test)

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

# Mostrar el reporte de clasificación del modelo SVM
print(classification_report(y_test, y_pred))

# Calcular el valor del área bajo la curva del gráfico ROC (AUC)
svc = RocCurveDisplay.from_estimator(rf, x_test, y_test, pos_label = "Alta")
plt.show()

# Guardar modelo entrenado
filename = "E:/ME/ALL/Project_NN_proteins/SCRIPTS/PYTHON/MODELS/idr_rf_model.sav"
joblib.dump(rf, filename)