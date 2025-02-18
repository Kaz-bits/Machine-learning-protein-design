# Cargar datos de la biblioteca (limpios)
df_library <- read.csv(file = "D:/FRET Script Biblioteca/data/FRET/all_FRET_biosensors.csv", 
                        header = TRUE, sep = ",")

# Cargar archivo de Sparrow
df_sparrow <- read.csv(file = "D:/ME/ALL/Project_NN_proteins/DATA/DATABASES/IDRBS_library_sparrow.csv", 
                       header = TRUE)

# Remover datos de SED1 de df_library
df_library <- df_library[!df_library$Construct == 201, ]


# Crear lista de constructos por remover
temp <- c(7, 16, 35, 76, 80, 87, 91, 92, 93, 98, 100, 103)

# Remover constructos
for (i in temp) {
  
  # Renovar data frame
  df_library <- df_library[!df_library$Construct== i, ]
  
}

# Obtener los nombres de los constructos
temp_names <- unique(df_library$Construct)

# Generar datos de delta FRET
temp_mean <- c()
temp_sd <- c()

for (a in temp_names) {
  
  # Obtener datos a 0 mM
  i <- df_library[df_library$Construct == a & 
                     df_library$Treatment == 0 & 
                     df_library$Plate == 1, ]$Normalized
  
  # Obtener datos a 1500 mM
  h <- df_library[df_library$Construct == a & 
                     df_library$Treatment == 1000, ]$Normalized
  
  # Obtener la media 
  temp_mean[length(temp_mean) + 1] <- mean(h - i)
  
  # Obtener la desviación estándar 
  temp_sd[length(temp_sd) + 1] <- sd(h - i)
  
}

# Constuir nuevo data frame
temp_delta <- data.frame("Construct" = temp_names,
                         "delta_fret" = temp_mean, 
                         "sd" = temp_sd)


# Asignar número de constructo a df_sparrow
df_sparrow$X <- seq(1:200)

# Filtrar datos de df_library con df_sparrow
df_sparrow <- df_sparrow[df_sparrow$X %in% df_library$Construct, ]

# Ordenar data frames anterior
df_sparrow <- df_sparrow[order(df_sparrow$X), ]
df_library <- df_library[order(df_library$Construct), ]

# Añadir a df_sparrow los datos de mean_delta de df_library
df_sparrow$delta <- temp_delta$delta_fret


# Obtener la cantidad de valores numéricos
temp <- nchar(temp_delta$Construct)

# Modificar nombres agregando "ceros"
temp_names <- ifelse(temp == 1, paste0("00", temp_delta$Construct), 
                     ifelse(temp == 2, paste0("0", temp_delta$Construct), 
                            temp_delta$Construct))

# Añadir nombres al data frame
df_sparrow$Construct <- temp_names

# Eliminar primer columna
df_sparrow$X <- NULL

# Agregar última columna a la primera
df_sparrow <- df_sparrow[, c(ncol(df_sparrow), 16, 1:15)]


# Guardar archivo
write.csv(x = df_sparrow, file = "D:/ME/ALL/Project_NN_proteins/DATA/DATABASES/IDRBS_library_sparrow_188.csv", 
          quote = FALSE, row.names = FALSE)


# Convertir valores negativos de delta FRET en cero
df_sparrow[df_sparrow$delta < 0, ]$delta <- 0

# Cambiar el valor de -1 de kappa por su valor real (calcular con CIDER)
df_sparrow[df_sparrow$kappa < 0, ]$kappa <- c(0.266, 0.233)



# Guardar archivo con datos limpiados
write.csv(x = df_sparrow, file = "D:/ME/ALL/Project_NN_proteins/DATA/DATABASES/IDRBS_library_sparrow_188_NN.csv", 
          quote = FALSE, row.names = FALSE)

