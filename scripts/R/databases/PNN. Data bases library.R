# Cargar archivo
fret <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_FRET_biosensors.csv", 
                 header = TRUE)

# Verificar la cantidad de biosensores
len <- length(names(table(fret$Construct)))

# Construir un data frame
temp_fret <- data.frame(matrix(nrow = len, ncol = 1))

# Cambiar el nombre de la primer columna
names(temp_fret)[1] <- "Construct"

# Agregar nombres de los constructos
temp_fret$Construct <- names(table(fret$Construct))

# Crear vector vacio
temp <- c() 


for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 1, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P1 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 2, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P2 <- temp

# Crear vector vacio
temp <- c() 

# Iterar sobre cada constructo
for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 200, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_200 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 400, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_400 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 600, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_600 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 800, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_800 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1000, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1000 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1500, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1500 <- temp


# Guardar archivo
write.table(x = temp_fret, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data_normalized.csv", 
            quote = FALSE, row.names = FALSE, sep = ",")







# Cargar archivo
fret <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_FRET_biosensors.csv", 
                 header = TRUE)

# Verificar la cantidad de biosensores
len <- length(names(table(fret$Construct)))

# Construir un data frame
temp_fret <- data.frame(matrix(nrow = len, ncol = 1))

# Cambiar el nombre de la primer columna
names(temp_fret)[1] <- "Construct"

# Agregar nombres de los constructos
temp_fret$Construct <- names(table(fret$Construct))

# Crear vector vacio
temp <- c() 


for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 1, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P1 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 2, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P2 <- temp

# Crear vector vacio
temp <- c() 

# Iterar sobre cada constructo
for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 200, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_200 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 400, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_400 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 600, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_600 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 800, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_800 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1000, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1000 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1500, ]$Normalized)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1500 <- temp


# Guardar archivo
write.table(x = temp_fret, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data_normalized.csv", 
            quote = FALSE, row.names = FALSE, sep = ",")








# Cargar archivo
fret <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_FRET_biosensors.csv", 
                 header = TRUE)

# Verificar la cantidad de biosensores
len <- length(names(table(fret$Construct)))

# Construir un data frame
temp_fret <- data.frame(matrix(nrow = len, ncol = 1))

# Cambiar el nombre de la primer columna
names(temp_fret)[1] <- "Construct"

# Agregar nombres de los constructos
temp_fret$Construct <- names(table(fret$Construct))

# Crear vector vacio
temp <- c() 


for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 1, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P1 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 2, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P2 <- temp

# Crear vector vacio
temp <- c() 

# Iterar sobre cada constructo
for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 200, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_200 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 400, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_400 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 600, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_600 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 800, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_800 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1000, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1000 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1500, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1500 <- temp


# Guardar archivo
write.table(x = temp_fret, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data_no_normalized.csv", 
            quote = FALSE, row.names = FALSE, sep = ",")







# Cargar archivo
fret <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_FRET_biosensors.csv", 
                 header = TRUE)

# Verificar la cantidad de biosensores
len <- length(names(table(fret$Construct)))

# Construir un data frame
temp_fret <- data.frame(matrix(nrow = len, ncol = 1))

# Cambiar el nombre de la primer columna
names(temp_fret)[1] <- "Construct"

# Agregar nombres de los constructos
temp_fret$Construct <- names(table(fret$Construct))

# Crear vector vacio
temp <- c() 


for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 1, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P1 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 0 & fret$Plate == 2, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_0_P2 <- temp

# Crear vector vacio
temp <- c() 

# Iterar sobre cada constructo
for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 200, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_200 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 400, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_400 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 600, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_600 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 800, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_800 <- temp


# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1000, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1000 <- temp

# Crear vector vacio
temp <- c() 

for (i in 1:nrow(temp_fret)) {
  
  # Extraer los datos por concentración y constructo 
  temp1 <- mean(fret[fret$Construct == i & fret$Treatment == 1500, ]$DxAm.DxDm)
  
  # Guardar datos de la media
  temp[length(temp) + 1] <- temp1  
  
}

# Añadir columna a temp_fret
temp_fret$NaCl_1500 <- temp


# Guardar archivo
write.table(x = temp_fret, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data_DxAm.DxDm.csv", 
            quote = FALSE, row.names = FALSE, sep = ",")


# Obtener datos entre bases de datos
temp_fret_1 <- temp_fret[as.numeric(temp_fret$Construct) %in% as.numeric(substr(df_sparrow$Entry, 7, 9)), ]

# Juntar data frames con df_sparrow y guardarlo en nueva variable
temp_fret_1 <- cbind(df_sparrow, temp_fret_1)

# Remover columnas y acomodar columnas
temp_fret_1 <- temp_fret_1[, -22]
temp_fret_1 <- temp_fret_1[, c(1:2, 20:21, 3:19, 22:29)]

# Cambiar nombre de columnas
names(temp_fret_1)[3] <- "delta_FRET"

# Guardar archivo
write.table(x = temp_fret_1, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data.csv", 
            quote = FALSE, row.names = FALSE, sep = ",")

