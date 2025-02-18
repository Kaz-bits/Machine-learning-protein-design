# Paqueterías
library(ggplot2)
library(dplyr)

# Cargar archivo
cider_parrot <- read.table(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/cider_library_PCA.csv", 
                           header = TRUE, sep = ",")

# Cargar biblioteca de biosensores
bios <- readxl::read_xlsx(path = "D:/ME/ALL/Project_IDP_D2P2/DATA/DATA FOR R/DATA FOR ANALYSIS/200 IDRBS Library.xlsx", 
                          sheet = 1, col_names = TRUE)


# Extraer los datos del cluster 1 y 2
cider_parrot_1 <- cider_parrot[cider_parrot$clustering == 1, ]
cider_parrot_2 <- cider_parrot[cider_parrot$clustering == 2, ]

# Añadir secuencias a los data frames de parrot
temp <- bios[as.numeric(substr(bios$Entry, start = 7, 9)) %in% cider_parrot_1$construct,]
temp1 <- bios[as.numeric(substr(bios$Entry, start = 7, 9)) %in% cider_parrot_2$construct,]

# Ordenar datos por constructo
temp <- temp[order(temp$Entry), ]
temp1 <- temp1[order(temp1$Entry), ]

# Ordenar datos de cider por constructo
cider_parrot_1 <- cider_parrot_1[order(cider_parrot_1$construct), ]
cider_parrot_2 <- cider_parrot_2[order(cider_parrot_2$construct), ]

# Añadir secuencias a cider_parrot_1
cider_parrot_1$sequence <- temp$`IDR sequence to order`
cider_parrot_2$sequence <- temp1$`IDR sequence to order`

# Ordenar datos
cider_parrot_1 <- cider_parrot_1[, c(1, 11, 2, 3:10)]
cider_parrot_2 <- cider_parrot_2[, c(1, 11, 2, 3:10)]

# Cambiar etiquetas de respuesta por numeros (alta = 2, media  = 1 y baja  = 0)
cider_parrot_1[cider_parrot_1$Response == "Alta", ]$Response <- 2
cider_parrot_1[cider_parrot_1$Response == "Media", ]$Response <- 1
cider_parrot_1[cider_parrot_1$Response == "Baja", ]$Response <- 0

cider_parrot_2[cider_parrot_2$Response == "Alta", ]$Response <- 2
cider_parrot_2[cider_parrot_2$Response == "Media", ]$Response <- 1
cider_parrot_2[cider_parrot_2$Response == "Baja", ]$Response <- 0

# Ordenar columnas
cider_parrot_1 <- cider_parrot_1[, c(1, 2, 10, 3:9, 11)]
cider_parrot_2 <- cider_parrot_2[, c(1, 2, 10, 3:9, 11)]

# Remover columnas
cider_parrot_1 <- cider_parrot_1[-c(4:11)]
cider_parrot_2 <- cider_parrot_2[-c(4:11)]

# Convertir datos de respuesta a enteros
cider_parrot_1$Response <- as.integer(cider_parrot_1$Response)
cider_parrot_2$Response <- as.integer(cider_parrot_2$Response)

# Guardar datos
write.table(x = cider_parrot_1, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/PARROT/cider_parrot_1.tsv", 
            quote = FALSE, row.names = FALSE, sep = " ")

write.table(x = cider_parrot_2, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/PARROT/cider_parrot_2.tsv", 
            quote = FALSE, row.names = FALSE, sep = " ")


# Observar distribución de los datos por respuesta del cluster 1
ggplot() +
  geom_histogram(data = cider_parrot_1, aes(x = mean_delta), 
                 bins = 30, fill = "steelblue", color = "black", size = 0.8) +
  theme_classic() +
  labs(x = expression(Delta*"FRET"), y = "Cantidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(-0.5, 1.8, 0.1)) +
  scale_y_continuous(expand = c(0, 0))


# Observar distribución de los datos por respuesta del cluster 2
ggplot() +
  geom_histogram(data = cider_parrot_2, aes(x = mean_delta), 
                 bins = 30, fill = "steelblue", color = "black", size = 0.8) +
  theme_classic() +
  labs(x = expression(Delta*"FRET"), y = "Cantidad") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(-0.5, 0.8, 0.1)) +
  scale_y_continuous(expand = c(0, 0))





# Cargar archivo
cider_parrot <- read.table(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/cider_library_PCA.csv", 
                           header = TRUE, sep = ",")

# Cargar biblioteca de biosensores
bios <- readxl::read_xlsx(path = "D:/ME/ALL/Project_IDP_D2P2/DATA/DATA FOR R/DATA FOR ANALYSIS/200 IDRBS Library.xlsx", 
                          sheet = 1, col_names = TRUE)

# Asegurar reproducibilidad
set.seed(1)

# Extraer una muestra de 149 datos
cider_parrot_sam <- sample_n(cider_parrot, size = 149)

# Verificar distribución de datos por respuesta
table(cider_parrot_sam$Response)

# Asegurar reproducibilidad
set.seed(2)

# Extraer 38 datos de media de forma aleatoria
temp_cider <- sample_n(cider_parrot_sam[cider_parrot_sam$Response == "Media", ], 35)

# Extraer 38 datos de alta de forma aleatoria
temp_cider1 <- sample_n(cider_parrot_sam[cider_parrot_sam$Response == "Alta", ], 35)

# Extraer los datos de baja
temp_cider2 <- cider_parrot_sam[cider_parrot_sam$Response == "Baja", ]

# Juntar data frames
temp_cider <- rbind(temp_cider, temp_cider1, temp_cider2)

# Añadir secuencias a los data frames de parrot
temp <- bios[as.numeric(substr(bios$Entry, start = 7, 9)) %in% temp_cider$construct,]

# Ordenar datos por constructo
temp <- temp[order(temp$Entry), ]

# Ordenar datos de cider por constructo
temp_cider <- temp_cider[order(temp_cider$construct), ]

# Añadir secuencias a cider_parrot_1
temp_cider$sequence <- temp$`IDR sequence to order`

# Ordenar datos
cider_parrot_sam <- temp_cider[, c(1, 11, 2, 3:10)]

# Cambiar etiquetas de respuesta por numeros (alta = 2, media  = 1 y baja  = 0)
cider_parrot_sam[cider_parrot_sam$Response == "Alta", ]$Response <- 2
cider_parrot_sam[cider_parrot_sam$Response == "Media", ]$Response <- 1
cider_parrot_sam[cider_parrot_sam$Response == "Baja", ]$Response <- 0

# Ordenar columnas
cider_parrot_sam <- cider_parrot_sam[, c(1, 2, 10, 3:9, 11)]

# Remover columnas
cider_parrot_sam <- cider_parrot_sam[-c(4:11)]

# Convertir datos de respuesta a enteros
cider_parrot_sam$Response <- as.integer(cider_parrot_sam$Response)

# Convertir identificadores a caracteres
cider_parrot_sam$construct <- as.character(cider_parrot_sam$construct)

# Guardar datos
write.table(x = cider_parrot_sam, file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/PARROT/cider_parrot_class.tsv", 
            quote = FALSE, row.names = FALSE, sep = " ")


table(cider_parrot_sam$Response)

