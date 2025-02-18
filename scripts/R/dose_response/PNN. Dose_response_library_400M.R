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

# Construir dos data frame vacíos
temp_df1 <- data.frame(matrix(ncol = 3))
temp_df2 <- data.frame(matrix(ncol = 1))

# Colocar nombres al data frame temporal
names(temp_df1)[1:3] <- c("NaCl_400", "Dose", "Construct")
names(temp_df2)[1] <- c("DxAm.DxDm")

for (a in list_names) {
  
  # Obtener data frame individual
  temp <- df_fret_200[df_fret_200$Construct == a, ]
  
  # Cambiar ceros por 0.01 para transformación logarítmica
  temp <- temp %>% mutate(NaCl_400 = ifelse((NaCl_400 == 0),
                                            yes = 0.01, 
                                            no = NaCl_400))
  
  # Calcular el promedio del valor de normalización para cada
  # tratamiento por cada réplica
  temp_summary <- temp %>% group_by(NaCl_400, Replicate, Construct) %>% 
    summarise("Mean" = mean(Normalized))
  
  # Predecir modelo no lineal para la concentración
  fit <- drm(formula = Mean ~ NaCl_400, 
             data = temp_summary, 
             fct = LL2.4())
  
  # Calcular la "dosis" de las IDRs al 10%, 50% y 90%
  dose <- round(ED(fit, c(10, 50, 90), interval = "fls")[1:3], digits = 2)
  
  # Generar data frame
  temp_df <- data.frame("NaCl_400" = dose, 
                        "Dose" = c("10", "50", "90"), 
                        "Construct" = a)
  
  # Unir data frames
  temp_df1 <- rbind(temp_df1, temp_df)
  
  
  # Predecir modelo no lineal para la la respuesta de FRET
  fit <- drm(formula = NaCl_400 ~ Mean, 
             data = temp_summary, 
             fct = LL2.4())
  
  # Calcular la "dosis" de las IDRs al 10%, 50% y 90%
  dose <- round(ED(fit, c(10, 50, 90), interval = "fls")[1:3], digits = 2)
  
  # Generar data frame
  temp_df3 <- data.frame("DxAm.DxDm" = dose)
  
  # Unir data frames
  temp_df2 <- rbind(temp_df2, temp_df3)
  
}

# Eliminar primer renglón
temp_df1 <- temp_df1[-1, ]
temp_df2 <- temp_df2[-1, ]

# Unir ambos data frames
df_dose <- cbind(temp_df1, temp_df2)

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
df_dose <- df_dose[, -c(5)]

# Reacomodar columnas
df_dose <- df_dose[, c(3, 2, 1, 4:33)]

# Renombrar columnas
names(df_dose)[4] <- "DxAm.DxDm"


# Generar data frames individuales (10, 50 y 90)
df_10 <- df_dose[df_dose$Dose == "10", ]
df_50 <- df_dose[df_dose$Dose == "50", ]
df_90 <- df_dose[df_dose$Dose == "90", ]


# Remover datos mayores a 10^3
df_10 <- df_10[df_10$NaCl_400 <= 507.28, ]
df_50 <- df_50[df_50$NaCl_400 <= 1458.58, ]
df_90 <- df_90[df_90$NaCl_400 <= 1473.97, ]

# Verificar los outliers de los datasets
grubbs.test(sort(df_50$NaCl_400), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$NaCl_400 == 0.05, ]
df_50 <- df_50[!df_50$NaCl_400 == 1458.58, ]

# Verificar los outliers de los datasets
grubbs.test(sort(df_50$NaCl_400), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$NaCl_400 == 4.86, ]
df_50 <- df_50[!df_50$NaCl_400 == 1247.96, ]


# Verificar los outliers de los datasets
grubbs.test(sort(df_50$NaCl_400), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$NaCl_400 == 8.57, ]
df_50 <- df_50[!df_50$NaCl_400 == 936.88, ]


# Verificar los outliers de los datasets
grubbs.test(sort(df_50$NaCl_400), type = 11)

# Remover outliers
df_50 <- df_50[!df_50$NaCl_400 == 75.87, ]
df_50 <- df_50[!df_50$NaCl_400 == 901.56, ]



# Gráficos de correlación [NaCl] (400 M) ----
# [NaCl] (400 M) vs kappa

# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$kappa, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$kappa, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = kappa, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(kappa), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.09, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.24, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_kappa.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_kappa.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# [NaCl] (400 M) vs kappa

# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$kappa1, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$kappa1, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = kappa1, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(kappa), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.08, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.20, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_kappa1.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_kappa1.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)



# [NaCl] (400 M) vs FCR
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$FCR, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$FCR, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = FCR, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "FCR", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.11, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.25, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_FCR.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_FCR.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs NCPR
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$NCPR, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$NCPR, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = NCPR, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "NCPR", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = -0.35, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = -0.21, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_NCPR.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_NCPR.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)






# [NaCl] (400 M) vs Hidropatía
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$hydropathy, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$hydropathy, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = hydropathy, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "Hidropatía", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 2.00, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 2.55, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_hidropatia.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_hidropatia.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs SCD
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$SCD, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$SCD, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = SCD, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "SCD", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = -0.5, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 10, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_SCD.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_SCD.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)




# [NaCl] (400 M) vs SHD
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$SHD, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$SHD, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = SHD, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "SHD", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 2.75, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 3.55, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_SHD.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_SHD.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs radio de giro
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$radius_of_gyration, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$radius_of_gyration, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = radius_of_gyration, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[g]), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 21, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 29, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_rg.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_rg.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs SHD
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$delta_FRET, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$delta_FRET, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = delta_FRET, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = expression(Delta*"FRET"), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.05, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.48, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_delta_FRET.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_delta_FRET.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs alifáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$aliphatic, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$aliphatic, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = aliphatic, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \nalifáticos", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.065, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.15, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_aliphatic.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_aliphatic.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs polares
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$polar, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$polar, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = polar, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \npolares", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.09, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.21, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_polar.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_polar.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs aromáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$aromatic, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$aromatic, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = aromatic, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = "Fracción de residuos \naromáticos", y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.006, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.029, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_aromatic.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_aromatic.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs aomáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$end_to_end_distance, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$end_to_end_distance, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = end_to_end_distance, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[e]), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 48, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 66, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_re.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_re.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)





# [NaCl] (400 M) vs aomáticos
# Valor de correlación
a <- round(unname(cor.test(df_50$NaCl_400, df_50$asphericity, method = "pearson")$estimate), 2)

# Valor de p
b <- format(cor.test(df_50$NaCl_400, df_50$asphericity, method = "pearson")$p.value, digits = 3)


# Gráfico 1. kappa vs FRET
ggplot(data = df_50, aes(x = asphericity, y = NaCl_400)) +
  geom_point() +
  theme_bw() +
  labs(x = bquote(R[e]), y = bquote(DxAm/DxDm[" [NaCl] (400M)"])) +
  geom_smooth(method = "lm") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank()) +
  annotate(geom = "text", y = 2.3, x = 0.355, label = paste("r = ", a)) +
  annotate(geom = "text", y = 2.3, x = 0.385, label = paste("p = ", b))

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_asph.pdf",
       device = "pdf", width = 4.5, height = 4, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "D:/ME/ALL/Project_NN_proteins/PLOTS/NaCl/400M/NaCl_400M_asph.png",
       device = "png", width = 4.5, height = 4, units = "in", dpi = 400)
