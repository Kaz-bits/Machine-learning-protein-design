# Cargar sparrow
from sparrow import read_fasta
import pandas as pd

# Cargar archivo fasta
protein_dictionary = read_fasta('my_fasta_file.fasta')

# Extraer nombres de las secuencias
temp_names = list(protein_dictionary)

# Crear una lista vacía
list_of_IDRs = []

# Obtener secuencias de las IDRs
protein_list = list(protein_dictionary.values())

# Generar lista vacía
temp_kappa = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_kappa.append(protein_list[i].kappa)


# Generar lista vacía
temp_FCR = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_FCR.append(protein_list[i].FCR)


# Generar lista vacía
temp_NCPR = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_NCPR.append(protein_list[i].NCPR)


# Generar lista vacía
temp_hydro = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_hydro.append(protein_list[i].hydrophobicity)


# Generar lista vacía
temp_SCD = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_SCD.append(protein_list[i].SCD)


# Generar lista vacía
temp_SHD = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_SHD.append(protein_list[i].SHD)


# Generar lista vacía
temp_aro = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_aro.append(protein_list[i].fraction_aromatic)
    

# Generar lista vacía
temp_ali = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_ali.append(protein_list[i].fraction_aliphatic)


# Generar lista vacía
temp_polar = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_polar.append(protein_list[i].fraction_polar)


# Generar lista vacía
temp_proline = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_proline.append(protein_list[i].fraction_proline)


# Generar lista vacía
temp_MW = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_MW.append(protein_list[i].molecular_weight)


# Generar lista vacía
temp_complex = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_complex.append(protein_list[i].complexity)


# Generar lista vacía
temp_rg = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_rg.append(protein_list[i].predictor.radius_of_gyration(use_scaled=True))


# Generar lista vacía
temp_re = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_re.append(protein_list[i].predictor.end_to_end_distance())


# Generar lista vacía
temp_prefactor = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_prefactor.append(protein_list[i].predictor.prefactor())



# Generar lista vacía
temp_scaling = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_scaling.append(protein_list[i].predictor.scaling_exponent())



# Generar lista vacía
temp_aes = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_aes.append(protein_list[i].predictor.asphericity())



# Generar lista vacía
temp_seq = []

# Iterar sobre cada secuencia
for i in range(len(protein_list)):
    temp_seq.append(protein_list[i].sequence)


# Crear una lista de listas con las variables generadas
sparrow_lists = temp_names, temp_seq, temp_kappa, temp_FCR, temp_rg, temp_NCPR, temp_hydro, temp_SCD, temp_SHD, temp_aro, temp_ali, temp_polar, temp_proline, temp_complex, temp_MW, temp_re, temp_prefactor, temp_scaling, temp_aes

# Crear el dataframe con el modulo de Pandas y transponer las columnas
# usando la función "T"
df_sparrow = pd.DataFrame(sparrow_lists, index = ["Construct", "sequence", "kappa", "FCR", "rg", "NCPR", "hydropathy", "SCD", "SHD", "aromatic", "aliphatic", "polar", "proline", "complex", "MW", "re", "prefactor", "scaling_exponent", "asphericity"]).T

# Guardar archivo en el directorio de elección
df_sparrow.to_csv("IDRBS_library_sparrow_200.csv")