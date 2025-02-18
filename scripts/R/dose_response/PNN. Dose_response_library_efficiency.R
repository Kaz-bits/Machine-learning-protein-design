# Paqueterías
library(ggplot2)
library(dplyr)
library(drc)

# Cargar archivo
temp <- read.csv(file = "C:/Users/Cesar/OneDrive/Escritorio/all_FRET_efficiency.csv", 
                 header = TRUE)


# Remover segundo cero
temp <- temp[!(temp$Plate == 2 & temp$Treatment == 0), ]

# Renumerar renglones
row.names(temp) <- seq(1, nrow(temp))

# Extraer todos los ceros
temp_zero <- temp[temp$Treatment == 0, ]

# Generar vector vacío
temp_mean <- c()

# Obtener todas las medias de cada constructo
for (construct in unique(temp$Construct)) {
  
  # Extraer datos por constructo
  temp_zero1 <- temp_zero[temp_zero$Construct == construct, ]
  
  # Extraer ceros por cada réplica
  i <- mean(temp_zero1[1:3, ]$Percentaje.efficiency)
  h <- mean(temp_zero1[4:6, ]$Percentaje.efficiency)
  j <- mean(temp_zero1[7:9, ]$Percentaje.efficiency)
  
  # Juntar medias en un vector
  temp_mean[length(temp_mean) + 1] <- i
  temp_mean[length(temp_mean) + 1] <- h
  temp_mean[length(temp_mean) + 1] <- j
  
}

# Generar repeticiones de las medias y añadirlo al data frame
temp$Mean <- rep(temp_mean, each = 21)

# Obtener el valor de eficiencia FRET normalizado
temp$Normalized <- (temp$FRET.FRET_D / temp$Mean)

# Guardar archivo
write.csv(x = temp, file = "D:/ME/ALL/Project_NN_proteins/DATA/DATABASES/all_FRET_efficiency.csv", 
          quote = FALSE, row.names = FALSE)

# Generar gráficos para cada constructo
for (i in unique(temp$Construct)) {
  
  # Construir gráficos de eficiencia FRET indivuales
  ggplot(data = temp[temp$Construct == i, ], 
         aes(x = Treatment/1000, y = Normalized)) +
    geom_point() +
    geom_smooth() +
    theme_classic() +
    labs(x = "[NaCl] (M)", y = "Normalized \n FRET efficiency") +
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12))
  
  #Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/PDF", paste0("FRET_efficiency_", i, ".pdf")),
         device = "pdf", width = 6, height = 5, units = "in", 
         dpi = 450) 
  
  #Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/PNG/", paste0("FRET_efficiency_", i, ".png")),
         device = "png", width = 6, height = 5, units = "in", 
         dpi = 450) 
  
}


# Gráfico dosis respuesta de las eficiencias FRET ----
# Construir modelo no lineal

# Obtener los nombres de los Constructos del
# archivo cargado
list_names <- names(table(temp$Construct))

# Realizar una iteración para obtener los datos de
# cada Constructo en la variable "temp"
temp_newdata <- data.frame()

for (a in list_names) {
  
  # Obtener data frame individual
  temp1 <- temp[temp$Construct == a, ]
  
  # Cambiar ceros por 0.01 para transformación logarítmica
  temp2 <- temp1 %>% mutate(Treatment = ifelse((Treatment == 0),
                                             yes = 0.01, 
                                             no = Treatment))
  
  # Calcular el promedio del valor de normalización para cada
  # tratamiento por cada réplica
  temp_summary <- temp2 %>% group_by(Treatment, Replicate, Construct) %>% 
    summarise("Mean" = mean(Normalized))
  
  # Predecir modelo no lineal
  fit <- drm(formula = Mean ~ Treatment, 
             data = temp_summary, 
             fct = LL2.4())
  
  # Crear nuevos datos para la curva del modelo
  newdata <- expand.grid(Treatment = exp(seq(0, 10, length = 100)))
  
  # Predecir intervalos de confianza
  pm <- predict(fit, newdata = newdata, interval = "confidence")
  
  # Datos nuevos con las predicciones al dataframe
  newdata$p <- pm[, 1]
  newdata$pmin <- pm[, 2]
  newdata$pmax <- pm[, 3]
  newdata$Construct <- a
  
  # Juntar todos los datos en un dataframe
  temp_newdata <- rbind(newdata, temp_newdata)
  
}



# Graficar la curva con los datos (colores)
ggplot() +
  geom_line(data = temp_newdata, 
            aes(x = (Treatment/1000), y = p, 
                group = Construct, color = Construct),
            size = 0.7, 
            lty = 1, alpha = 0.3, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nFRET efficiency") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(0.95, 3)) +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25), expand = c(0,0))


#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/DRC_all_bios_efficiency.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/DRC_all_bios_efficiency.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 



# Graficar la curva con los datos de 163 en rojo y SED1 en azul
ggplot() +
  geom_line(data = temp_newdata, 
            aes(x = (Treatment/1000), y = p, 
                group = Construct, color = Construct),
            size = 0.7, color = "gray", 
            lty = 1, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = temp_newdata[temp_newdata$Construct == 163, ], 
            aes(x = (Treatment/1000), y = p, 
                group = Construct),
            size = 1.2, color = "red", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  geom_line(data = temp_newdata[temp_newdata$Construct == 201, ], 
            aes(x = (Treatment/1000), y = p, 
                group = Construct),
            size = 1.2, color = "blue", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nFRET efficiency") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(0.95, 3)) +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25),
                     expand = c(0,0))

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/DRC_all_bios_efficiency_CAHS_SED1.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/FRET EFFICIENCY/DRC_all_bios_efficiency_CAHS_SED1.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 




