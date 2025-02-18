# Paqueterías
library(ggplot2)
library(forcats)

# Load data
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", 
                       header = TRUE)


# Specifiy path for individual NCPR files 
temp_path <- "D:/MASTER_FILES/DATA/NCPR_Constanza"
  
# List of biosensors
temp_list <- df_sparrow$Construct

# Extract high response biosensors
temp_list_high <- paste0(df_sparrow[df_sparrow$Response == "High", ]$Construct, ".csv")

# Extract medium response biosensors
temp_list_med <- paste0(df_sparrow[df_sparrow$Response == "Medium", ]$Construct, ".csv")

# Extract low response biosensors
temp_list_low <- paste0(df_sparrow[df_sparrow$Response == "Low", ]$Construct, ".csv")

# Data frame vacío
ncpr <- data.frame(matrix(nrow = 1, ncol = 6))

# Juntar los datos de biosensores de respuesta alta
for (i in temp_list_high) {
  
  # Cargar archivo
  temp <- read.csv(file = file.path(temp_path, i))
  
  # Colocar nombres en ncpr
  names(ncpr) <- names(temp)
  
  # Juntar datos
  ncpr <- rbind(ncpr, temp)
  
}

# Eliminar primer renglón
ncpr <- ncpr[-1, ]

# Calcular la cantidad de residuos positivos
pos <- length(ncpr[ncpr$Charge == "Positivo", ]$Charge)

# Calcular la cantidad de residuos negativos
neg <- length(ncpr[ncpr$Charge == "Negativo", ]$Charge)

# Calcular la cantidad de residuos no cargados
noch <- length(ncpr[ncpr$Charge == "No charge", ]$Charge)

# Juntar vectores
a <- c(pos, neg, noch)

# Gráfico con cantida de cargas
high <- ggplot(data = ncpr) +
  geom_bar(aes(x = Charge, y = ..count../sum(..count..), 
               fill = Charge), 
           color = "black", size = 0.8, width = 0.5,
           show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Frequency") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#FF6464", "#CFCFCF", "#3771FF"))


high

# Guardar gráfico
ggsave(plot = high, filename = "D:/MASTER_FILES/PLOTS/charge_dist_high.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = high, filename = "D:/MASTER_FILES/PLOTS/charge_dist_high.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Data frame vacío
ncpr <- data.frame(matrix(nrow = 1, ncol = 6))

# Juntar los datos de biosensores de respuesta alta
for (i in temp_list_med) {
  
  # Cargar archivo
  temp <- read.csv(file = file.path(temp_path, i))
  
  # Colocar nombres en ncpr
  names(ncpr) <- names(temp)
  
  # Juntar datos
  ncpr <- rbind(ncpr, temp)
  
}

# Eliminar primer renglón
ncpr <- ncpr[-1, ]

# Calcular la cantidad de residuos positivos
pos <- length(ncpr[ncpr$Charge == "Positivo", ]$Charge)

# Calcular la cantidad de residuos negativos
neg <- length(ncpr[ncpr$Charge == "Negativo", ]$Charge)

# Calcular la cantidad de residuos no cargados
noch <- length(ncpr[ncpr$Charge == "No charge", ]$Charge)

# Juntar vectores
a <- c(pos, neg, noch)

# Gráfico con cantida de cargas
med <- ggplot(data = ncpr) +
  geom_bar(aes(x = Charge, y = ..count../sum(..count..),
               fill = Charge), 
           color = "black", size = 0.8, width = 0.5,
           show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Frequency") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#FF6464", "#CFCFCF", "#3771FF"))

med

# Guardar gráfico
ggsave(plot = med, filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/obs/PDF/charge_dist_medium.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = med, filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/obs/PNG/charge_dist_medium.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)




# Data frame vacío
ncpr <- data.frame(matrix(nrow = 1, ncol = 6))

# Juntar los datos de biosensores de respuesta alta
for (i in temp_list_low) {
  
  # Cargar archivo
  temp <- read.csv(file = file.path(temp_path, i))
  
  # Colocar nombres en ncpr
  names(ncpr) <- names(temp)
  
  # Juntar datos
  ncpr <- rbind(ncpr, temp)
  
}

# Eliminar primer renglón
ncpr <- ncpr[-1, ]

# Calcular la cantidad de residuos positivos
pos <- length(ncpr[ncpr$Charge == "Positivo", ]$Charge)

# Calcular la cantidad de residuos negativos
neg <- length(ncpr[ncpr$Charge == "Negativo", ]$Charge)

# Calcular la cantidad de residuos no cargados
noch <- length(ncpr[ncpr$Charge == "No charge", ]$Charge)

# Juntar vectores
a <- c(pos, neg, noch)

# Gráfico con cantida de cargas
low <- ggplot(data = ncpr) +
  geom_bar(aes(x = Charge, y = ..count../sum(..count..), 
               fill = Charge), 
           color = "black", size = 0.8, width = 0.5,
           show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Frequency") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#FF6464", "#CFCFCF", "#3771FF"))

low

# Guardar gráfico
ggsave(plot = low, filename = "D:/MASTER_FILES/PLOTS/charge_dist_low.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(plot = low, filename = "D:/MASTER_FILES/PLOTS/charge_dist_low.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)

