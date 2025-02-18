# Paqueterías ----
library(ggplot2)
library(dplyr)
library(drc)
library(broom)
library(scales)
library(outliers)

# Cargar archivo ----
df_fret_200 <- read.csv(file = "D:/FRET Script Biblioteca/data/FRET/all_FRET_biosensors.csv", 
                        header = TRUE, sep = ",")

# Cargar datos de cider de la biblioteca
df_cider <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/all_fret_data.csv", 
                     header = TRUE)

# Cargar datos de cider previos de la biblioteca
df_cider1 <- read.csv(file = "D:/ME/ALL/Project_IDP_D2P2/DATA/BASES DE DATOS/df_fret_hml_delta.csv", 
                     header = TRUE)

# Ordenar datos de df_cider y df_cider 1
df_cider <- df_cider[order(df_cider$Entry), ]
df_cider1 <- df_cider1[order(df_cider1$construct), ]

# Añadir datos de kappa de df_cider1 en df_cider
df_cider$kappa1 <- df_cider1$kappa

# Remover los datos de kappa de -1
df_cider <- df_cider[!(df_cider$kappa == -1), ]



# Construir modelo no lineal ----

# Obtener los nombres de los Constructos del
# archivo cargado
list_names <- names(table(df_fret_200$Construct))

# Realizar una iteración para obtener los datos de
# cada Constructo en la variable "temp"
temp_newdata <- data.frame()

for (a in list_names) {
  
  # Obtener data frame individual
  temp <- df_fret_200[df_fret_200$Construct == a, ]
  
  # Cambiar ceros por 0.01 para transformación logarítmica
  temp <- temp %>% mutate(Treatment = ifelse((Treatment == 0),
                                             yes = 0.01, 
                                             no = Treatment))
  
  # Calcular el promedio del valor de normalización para cada
  # tratamiento por cada réplica
  temp_summary <- temp %>% group_by(Treatment, Replicate, Construct) %>% 
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

# Graficar la curva con los datos con SED1 en azul
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
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nDxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(0.95, 3)) +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25),
                     expand = c(0,0))

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/DRC_all_bios_CAHS.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/DRC_all_bios_CAHS.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 



# Graficar la curva con los datos (colores)
ggplot() +
  geom_line(data = temp_newdata, 
            aes(x = (Treatment/1000), y = p, 
                group = Construct, color = Construct),
            size = 0.7, 
            lty = 1, alpha = 0.3, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nDxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(0.95, 3)) +
  scale_x_continuous(expand = c(0,0))


#Guardar gráfico
ggsave(filename = "D:/FRET Script Biblioteca/data/FRET/DRC_all_bios_colores.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/FRET Script Biblioteca/data/FRET/DRC_all_bios_colores.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 



# Graficar la curva con los datos de CAHS y P53
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
  geom_line(data = temp_newdata[temp_newdata$Construct == 7, ], 
            aes(x = (Treatment/1000), y = p, 
                group = Construct),
            size = 1.2, color = "blue", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nDxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(0.95, 3)) +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25),
                     expand = c(0,0))

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/DRC_all_bios_p53_CAHS.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NEW_DATA/DRC_all_bios_p53_CAHS.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 



# Graficar factores de transcripción ----
# Lista de factores de transcripción
ft_list <- as.numeric(c("45", "54", "64", "77", "81", "82", "136", 
             "138", "139", "140", "141", "142", "143", "144", 
             "145", "146", "147", "148", "149", "165", "137"))

# Crear data frame vacío
temp_df1 <- data.frame(matrix(nrow = 1, ncol = ncol(temp_newdata)))

# Colocar nombres de las columnas
names(temp_df1) <- names(temp_newdata)

for (i in ft_list) {
  
  # Extraer datos de temp_newdata
  temp <- temp_newdata[temp_newdata$Construct == i, ]
  
  # Juntarlos en un data frame
  temp_df1 <- rbind(temp_df1, temp)

}

# Remover columna 1
temp_df1 <- temp_df1[-1, ]

# Gráfico de 21 factores de transcripción
ggplot() +
  geom_line(data = temp_df1, 
            aes(x = (Treatment/1000), y = p, 
                group = Construct, color = Construct), size = 0.7, 
            lty = 1, alpha = 0.8, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized \nDxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(1, 1.8)) +
  scale_x_continuous(expand = c(0, 0))

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_IDP_D2P2/PLOTS/FRET/Dose_response/DRC_all_factors_colors.pdf",
       device = "pdf", width = 6, height = 5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_IDP_D2P2/PLOTS/FRET/Dose_response/DRC_all_factors_colors.png",
       device = "png", width = 6, height = 5, units = "in", 
       dpi = 450) 


# Gráfico de factores de transcripción 144 (rojo) y 045 (azul)
ggplot() +
  geom_line(data = temp_df1, 
            aes(x = (Treatment/1000), y = p, 
                group = Construct, color = Construct),
            size = 0.7, color = "gray", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  geom_line(data = temp_df1[temp_df1$Construct == 144, ], 
            aes(x = (Treatment/1000), y = p, 
                group = Construct),
            size = 0.8, color = "red", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  geom_line(data = temp_df1[temp_df1$Construct == 45, ], 
            aes(x = (Treatment/1000), y = p, 
                group = Construct),
            size = 0.8, color = "blue", 
            lty = 1, alpha = 0.6, show.legend = FALSE) +
  theme_classic() +
  labs(x = "[NaCl] (M)", 
       y = "Normalized DxAm/DxDm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(xlim = c(0, 1.55), 
                  ylim = c(1, 1.8)) +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25),
                     expand = c(0, 0))


#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_IDP_D2P2/PLOTS/FRET/Dose_response/DRC_all_factors_144_045.pdf",
       device = "pdf", width = 4, height = 3.5, units = "in", 
       dpi = 450) 

#Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_IDP_D2P2/PLOTS/FRET/Dose_response/DRC_all_factors_144_045.png",
       device = "png", width = 4, height = 3.5, units = "in", 
       dpi = 450) 




# Construir modelo no lineal para ED50 ----
# Obtener los nombres de los Constructos del
# archivo cargado
list_names <- names(table(df_fret_200$Construct))

# Realizar una iteración para obtener los datos de
# cada Constructo en la variable "temp"
temp_newdata <- data.frame()

# Construir dos data frame vacíos
temp_df1 <- data.frame(matrix(ncol = 3))
temp_df2 <- data.frame(matrix(ncol = 1))

# Colocar nombres al data frame temporal
names(temp_df1)[1:3] <- c("Treatment", "Dose", "Construct")
names(temp_df2)[1] <- c("DxAm.DxDm")

# Crear vector vacio
temp_slope <- c()

for (a in list_names) {
  
  # Obtener data frame individual
  temp <- df_fret_200[df_fret_200$Construct == a, ]
  
  # Cambiar ceros por 0.01 para transformación logarítmica
  temp <- temp %>% mutate(Treatment = ifelse((Treatment == 0),
                                             yes = 0.01, 
                                             no = Treatment))
  
  # Calcular el promedio del valor de normalización para cada
  # tratamiento por cada réplica
  temp_summary <- temp %>% group_by(Treatment, Replicate, Construct) %>% 
    summarise("Mean" = mean(Normalized))
  
  # Predecir modelo no lineal para la concentración
  fit <- drm(formula = Mean ~ Treatment, 
             data = temp_summary, 
             fct = LL2.4(names = c("Hill slope", "Min", "Max", "EC50")))
  
  # Calcular la "dosis" de las IDRs al 10%, 50% y 90%
  dose <- round(ED(fit, c(10, 50, 90), interval = "fls")[1:3], digits = 2)
  
  # Generar data frame con las pendientes de cada constructo
  temp_slope[length(temp_slope) + 1] <- c(fit$fit$par[1])
  
  # Generar data frame
  temp_df <- data.frame("Treatment" = dose, 
                        "Dose" = c("10", "50", "90"), 
                        "Construct" = a)
  
  # Unir data frames
  temp_df1 <- rbind(temp_df1, temp_df)
  
  
  # Predecir modelo no lineal para la la respuesta de FRET
  fit <- drm(formula = Treatment ~ Mean, 
             data = temp_summary, 
             fct = LL2.4())
  
  # Calcular la "dosis" de las IDRs al 10%, 50% y 90%
  dose <- round(ED(fit, c(10, 50, 90), interval = "fls")[1:3], digits = 2)
  
  
  
  # Generar data frame
  temp_df3 <- data.frame("DxAm.DxDm" = dose)
  
  # Unir data frames
  temp_df2 <- rbind(temp_df2, temp_df3)
  
}

# Generar los valores de slope para todos los constructos
temp_slope <- rep(temp_slope, each = 3)

# Eliminar primer renglón
temp_df1 <- temp_df1[-1, ]
temp_df2 <- temp_df2[-1, ]

# Unir ambos data frames
df_dose <- cbind(temp_df1, temp_df2)

# Añadir columna de slopes a df_dose
df_dose$slope <- temp_slope

# Transformar los IDs de df_cider
df_cider$Entry <- as.numeric(substr(df_cider$Entry, 7, 9))

# Transformar IDs de df_dose
df_dose$Construct <- as.numeric(df_dose$Construct)

# Extraer datos de CIDER
df_dose <- df_dose[df_dose$Construct %in% df_cider$Entry, ]

# Guardar data frame de cider en nueva variable
temp <- df_cider

# Añadir temp a df_cider
df_cider <- rbind(df_cider, temp)

# Añadir nuevamente temp a df_cider
df_cider <- rbind(df_cider, temp)

# Ordenar datos de df_cider por constructo
df_cider <- df_cider[order(df_cider$Entry), ]

# Añadir df_cider a df_dose
df_dose <- cbind(df_dose, df_cider)

# Remover columnas
df_dose <- df_dose[, -c(6)]

# Reacomodar columnas
df_dose <- df_dose[, c(3, 2, 1, 4:34)]

# Renombrar columnas
names(df_dose)[4] <- "DxAm.DxDm"


# Generar data frames individuales (10, 50 y 90)
df_10 <- df_dose[df_dose$Dose == "10", ]
df_50 <- df_dose[df_dose$Dose == "50", ]
df_90 <- df_dose[df_dose$Dose == "90", ]


# Remover datos mayores a 10^3
df_10 <- df_10[df_10$Treatment <= 507.28, ]
df_50 <- df_50[df_50$Treatment <= 1458.58, ]
df_90 <- df_90[df_90$Treatment <= 1473.97, ]

# Verificar los outliers de los datasets
grubbs.test(sort(df_50$Treatment), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$Treatment == 0.05, ]
df_50 <- df_50[!df_50$Treatment == 1458.58, ]

# Verificar los outliers de los datasets
grubbs.test(sort(df_50$Treatment), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$Treatment == 4.86, ]
df_50 <- df_50[!df_50$Treatment == 1247.96, ]


# Verificar los outliers de los datasets
grubbs.test(sort(df_50$Treatment), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$Treatment == 8.57, ]
df_50 <- df_50[!df_50$Treatment == 936.88, ]


# Verificar los outliers de los datasets
grubbs.test(sort(df_50$Treatment), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$Treatment == 75.87, ]
df_50 <- df_50[!df_50$Treatment == 901.56, ]


# Gráficos de correlación ED50 ----
# ED50 vs kappa

# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$kappa, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$kappa, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = kappa, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(kappa), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.08, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.22, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_kappa.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_kappa.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# ED50 vs kappa

# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$kappa1, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$kappa1, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = kappa1, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(kappa), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.08, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.20, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_kappa1.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_kappa1.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# ED50 vs FCR
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$FCR, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$FCR, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = FCR, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "FCR", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.11, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.26, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_FCR.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)
 
# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_FCR.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs NCPR
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$NCPR, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$NCPR, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = NCPR, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "NCPR", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = -0.35, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = -0.20, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_NCPR.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_NCPR.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)






# ED50 vs Hidropatía
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$hydropathy, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$hydropathy, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = hydropathy, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "Hidropatía", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 2.08, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 2.70, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_hidropatia.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_hidropatia.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)







# ED50 vs Desorden promovido
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$Disorder.promoting, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$Disorder.promoting, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = Disorder.promoting, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "Desorden promovido", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.66, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.73, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_disorder.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_disorder.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# ED50 vs SCD
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$SCD, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$SCD, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = SCD, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "SCD", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = -0.5, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 9, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_SCD.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_SCD.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# ED50 vs SHD
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$SHD, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$SHD, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = SHD, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "SHD", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 2.75, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 3.55, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_SHD.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_SHD.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs radio de giro
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$radius_of_gyration, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$radius_of_gyration, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = radius_of_gyration, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[g]), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 21, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 28, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_rg.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_rg.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs SHD
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$delta_FRET, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$delta_FRET, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = delta_FRET, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(Delta*"FRET"), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.05, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.5, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_delta_FRET.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_delta_FRET.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs alifáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$aliphatic, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$aliphatic, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = aliphatic, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \nalifáticos", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.065, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.15, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_aliphatic.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_aliphatic.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs polares
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$polar, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$polar, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = polar, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \npolares", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.09, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.22, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_polar.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_polar.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs aromáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$aromatic, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$aromatic, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = aromatic, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \naromáticos", y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.006, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.028, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_aromatic.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_aromatic.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs aomáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$end_to_end_distance, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$end_to_end_distance, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = end_to_end_distance, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[e]), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 48, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 66, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_re.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_re.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# ED50 vs aomáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$asphericity, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$asphericity, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = asphericity, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[e]), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.355, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.385, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_asph.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_asph.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# ED50 vs aomáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$Treatment, df_50$NaCl_0_P1, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$Treatment, df_50$NaCl_0_P1, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = NaCl_0_P1, y = Treatment/1000)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(DxAm/DxDm["[NaCl] (0M)"]), y = "ED50") +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2, x = 0.65, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2, x = 0.95, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_NaCl_0M.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DOSE/ED50/ED50_NaCl_0M.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# Extraer datos de ED50 con todas las concentraciones
df_50[, c(26, 28:33)]

# Gráficos smooth ----
# Construir data frame vacío
temp_scaled <- data.frame(matrix(nrow = 1, ncol = ncol(df_fret_200) + 1))

# Agregar nombres
names(temp_scaled) <- c(names(df_fret_200), "scaled_y")

for (i in unique(df_fret_200$Construct)) {
  
  # Extraer datos del constructo
  temp <- df_fret_200[df_fret_200$Construct == i, ]
  
  # Escalar los datos de FRET de 0 a 100 para cada constructo
  scaled_y_matrix <- scale(temp$Normalized, 
        center = min(temp$Normalized), 
        scale = max(temp$Normalized) - min(temp$Normalized)) * 100
  
  # Extraer los datos escalados
  scaled_y <- scaled_y_matrix[1:nrow(scaled_y_matrix), 1]
  
  # Añadir columna
  temp$scaled_y <- scaled_y
  
  # Añadir datos al data frame
  temp_scaled <- rbind(temp_scaled, temp)
  
}


# Eliminar primer renglón
temp_scaled <- temp_scaled[-1, ]

# Calcular el promedio del valor de normalización para cada
# tratamiento por cada réplica
temp_summary <- temp_scaled %>% group_by(Treatment, Replicate, Construct) %>% 
  summarise("Mean" = mean(scaled_y))

# Generar gráficos con el eje Y estático
for (i in unique(temp_scaled$Construct)) {
  
  # Obtener los datos para cada constructo
  temp <- temp_summary[temp_summary$Construct == i, ]
  
  # Crear nombre de cada gráfico
  name <- paste0("IDRBS-", i)
  
  # Construir gráfico con smooth
  ggplot(data = temp, aes(x = Treatment / 1000, y = Mean)) + 
    geom_point() +
    geom_smooth() +
    theme_classic() +
    labs(x = "[NaCl] (mM)", y = "Respuesta de FRET") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    coord_cartesian(ylim = c(0, 120), xlim = c(0, 1.5)) +
    scale_y_continuous(breaks = seq(0, 120, 30)) +
    annotate(geom = "text", x = 0.75, y = 115, label = name)
  
  # Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL/PDF/", paste0(name, ".pdf")),
         device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)
  
  # Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL/PNG/", paste0(name, ".png")),
         device = "png", width = 4.5, height = 4, units = "in", dpi = 400)
  
  
  # Construir gráfico con smooth con línea al 50% de respuesta
  ggplot(data = temp, aes(x = Treatment / 1000, y = Mean)) + 
    geom_point() +
    geom_smooth() +
    geom_hline(yintercept = 77.715, lty = 2) +
    geom_hline(yintercept = 50, lty = 2) +
    geom_hline(yintercept = 17.041, lty = 2) +
    theme_classic() +
    labs(x = "[NaCl] (mM)", y = "Respuesta de FRET") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    coord_cartesian(ylim = c(0, 120), xlim = c(0, 1.5)) +
    scale_y_continuous(breaks = seq(0, 120, 30)) +
    annotate(geom = "text", x = 0.75, y = 115, label = name)
  
  # Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL/PDF/", paste0(name, ".pdf")),
         device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)
  
  # Guardar gráfico
  ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL/PNG/", paste0(name, ".png")),
         device = "png", width = 4.5, height = 4, units = "in", dpi = 400)
  
  
}




# Construir gráfico con smooth con línea al 50% de respuesta
ggplot(data = temp, aes(x = Treatment / 1000, y = Mean)) + 
  geom_point() +
  geom_smooth() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 17.041, ymax = 77.715), 
            alpha = 0.01, fill = "#A4D3EE") +
  theme_classic() +
  labs(x = "[NaCl] (mM)", y = "Respuesta de FRET") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 120), xlim = c(0, 1.5)) +
  scale_y_continuous(breaks = seq(0, 120, 30)) +
  annotate(geom = "text", x = 0.75, y = 115, label = name)

# Guardar gráfico
ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL_RANGE/PDF/", paste0(name, ".pdf")),
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = file.path("D:/ME/ALL/Project_NN_proteins/PLOTS/SMOOTH/INDIVIDUAL_RANGE/PNG/", paste0(name, ".png")),
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# Obtención del delta ED50 - FRET 0 mM ----
temp_change <- ((df_50$DxAm.DxDm - df_50$NaCl_0_P1) / (df_50$NaCl_0_P1)) * 100

# Crear vector vacío
temp_norm <- c()

# Normalizar los datos
for (i in temp_change) {
  
  # Normalizar
  temp_norm[length(temp_norm) + 1] <- ((i - min(temp_change)) / (max(temp_change) - min(temp_change)) * 100)

}

# Añadir vector temp_norm al data frame df_50
df_50$change <- temp_norm

# Ordenar los datos por respuesta de FRET
df_50$Response <- factor(df_50$Response, levels = c("Alta", "Media", "Baja"))

# Graficar la distribución de respuesta con ED50
ggplot(data = df_50, aes(x = Response, y = change, fill = Response)) +
  geom_boxplot(show.legend = FALSE, width = 0.5) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) + 
  labs(x = "Respuesta", y = "Sensibilidad de respuesta") + 
  scale_fill_manual(values = c("#3288FF", "#65BFFF", "#99E5FF"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/delta_ed50_change.pdf",
device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/delta_ed50_change.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)



# Calcular la magnitud del cambio
temp_change <- ((df_50$NaCl_1500 - df_50$NaCl_0_P1) / (df_50$NaCl_0_P1)) * 100

# Crear vector vacío
temp_norm <- c()

# Normalizar los datos de delta FRET
for (i in temp_change) {
  
  # Normalizar
  temp_norm[length(temp_norm) + 1] <- ((i - min(temp_change)) / (max(temp_change) - min(temp_change)) * 100)
  
}

# Añadir vector temp_norm al data frame df_50
df_50$change_delta <- temp_norm

# Graficar la distribución de respuesta con delta FRET
ggplot(data = df_50, aes(x = Response, y = change_delta, fill = Response)) +
  geom_boxplot(show.legend = FALSE, width = 0.5) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) + 
  labs(x = "Respuesta", y = "Estado conformacional") + 
  scale_fill_manual(values = c("#3288FF", "#65BFFF", "#99E5FF"))


# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/delta_fret_change.pdf",
       device = "pdf", width = 3.5, height = 3, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/delta_fret_change.png",
       device = "png", width = 3.5, height = 3, units = "in", dpi = 400)


# Determinar outliers de ED50 (change)
grubbs.test(sort(df_50$change), type = 11)

# Remover outliers de ED50 (change)
df_50 <- df_50[!df_50$change == 0, ]
df_50 <- df_50[!df_50$change == 100, ]

# Determinar outliers de ED50 (change)
grubbs.test(sort(df_50$change), type = 11)

# Remover outliers de ED50 (change)
df_50 <- df_50[!df_50$change == min(df_50$change), ]
df_50 <- df_50[!df_50$change == max(df_50$change), ]



# Determinar outliers de delta FRET (change_delta)
grubbs.test(sort(df_50$change_delta), type = 11)

# Remover outliers de ED50 (change)
df_50 <- df_50[!df_50$change_delta == min(df_50$change_delta), ]
df_50 <- df_50[!df_50$change_delta == max(df_50$change_delta), ]


# Correlación entre ED50 (change) vs Delta FRET ----
# Valor de correlación
a <- round(unname(cor.test(df_50$DxAm.DxDm, df_50$delta_FRET, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$DxAm.DxDm, df_50$delta_FRET, method = "pearson")$p.value, digits = 3)


# Gráfico de correlación entre ED50 (sensibilidad) vs delta FRET (estado conformacional)
ggplot(data = df_50, aes(x = delta_FRET, y = Treatment / 1000)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Compactación (ΔFRET)", y = "Sensibilidad (ED50)") +
  ggtitle(label = paste("r = ", a, "p = ", b)) +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        panel.grid = element_blank())

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/ed50_vs_delta.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/DELTA/ed50_vs_delta.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)


# Guardar archivo de df_50 (sin outliers)
write.csv(x = df_50, file = "D:/ME/ALL/Project_NN_proteins/DATA/DATABASES/ED50/ED50_all_biosensors_wth_outliers.csv", 
          quote = FALSE, row.names = FALSE)




