# Paqueterías
library(ggplot2)

# Cargar datos de sparrow
df_sparrow <- read.csv(file = "D:/MASTER_FILES/DATA/IDRBS_library_sparrow_188.csv", 
                       header = TRUE)

# Obtener las estadísticas sin los outliers
df_sparrow <- df_sparrow[!df_sparrow$delta < 0, ]

# Obtener el valor del tercer cuartil
a <- summary(df_sparrow$delta)[[5]]

# Obtener el valor del primer cuartil
b <- summary(df_sparrow$delta)[[2]]

# Obtener los datos de sparrow que tiene valores de delta FRET mayores a "a"
df_sparrow_a <- df_sparrow[df_sparrow$delta >= a, ]

# Obtener los datos de sparrow que tiene valores de delta FRET menores a "b"
df_sparrow_b <- df_sparrow[df_sparrow$delta < b, ]

# Obtener los datos de "df_sparrow_a" en "df_p_value"
df_p_value_a <- df_p_value[as.numeric(df_p_value$Construct) %in% df_sparrow_a$Construct, ]

# Obtener los datos de "df_sparrow_b" en "df_p_value"
df_p_value_b <- df_p_value[as.numeric(df_p_value$Construct) %in% df_sparrow_b$Construct, ]


# Transformar datos de hidropatía a una escala menor
df_sparrow_a$hydropathy <- (df_sparrow_a$hydropathy) / 10
df_sparrow_b$hydropathy <- (df_sparrow_b$hydropathy) / 10

# Transformar datos de Rg a una escala menor
df_sparrow_a$rg <- (df_sparrow_a$rg) / 100
df_sparrow_b$rg <- (df_sparrow_b$rg) / 100

# Transformar datos de SCD a una escala menor
df_sparrow_a$SCD <- (df_sparrow_a$SCD) / 10
df_sparrow_b$SCD <- (df_sparrow_b$SCD) / 10


# Construir data frames con los datos de sparrow_a
temp_df1a <- df_sparrow_a[, c(1, 7)]
temp_df2a <- df_sparrow_a[, c(1, 8)]
temp_df3a <- df_sparrow_a[, c(1, 9)]
temp_df4a <- df_sparrow_a[, c(1, 10)]
temp_df5a <- df_sparrow_a[, c(1, 19)]
temp_df6a <- df_sparrow_a[, c(1, 3)]
temp_df7a <- df_sparrow_a[, c(1, 14)]
temp_df8a <- df_sparrow_a[, c(1, 15)]
temp_df9a <- df_sparrow_a[, c(1, 16)]
temp_df10a <- df_sparrow_a[, c(1, 11)]


# Construir data frames con los datos de sparrow_b
temp_df1b <- df_sparrow_b[, c(1, 7)]
temp_df2b <- df_sparrow_b[, c(1, 8)]
temp_df3b <- df_sparrow_b[, c(1, 9)]
temp_df4b <- df_sparrow_b[, c(1, 10)]
temp_df5b <- df_sparrow_b[, c(1, 19)]
temp_df6b <- df_sparrow_b[, c(1, 3)]
temp_df7b <- df_sparrow_b[, c(1, 14)]
temp_df8b <- df_sparrow_b[, c(1, 15)]
temp_df9b <- df_sparrow_b[, c(1, 16)]
temp_df10b <- df_sparrow_b[, c(1, 11)]

# Cambiar nombres de las columnas por uno común
names(temp_df1a)[2] <- "Value"
names(temp_df2a)[2] <- "Value"
names(temp_df3a)[2] <- "Value"
names(temp_df4a)[2] <- "Value"
names(temp_df5a)[2] <- "Value"
names(temp_df6a)[2] <- "Value"
names(temp_df7a)[2] <- "Value"
names(temp_df8a)[2] <- "Value"
names(temp_df9a)[2] <- "Value"
names(temp_df10a)[2] <- "Value"

# Cambiar nombres de las columnas por uno común
names(temp_df1b)[2] <- "Value"
names(temp_df2b)[2] <- "Value"
names(temp_df3b)[2] <- "Value"
names(temp_df4b)[2] <- "Value"
names(temp_df5b)[2] <- "Value"
names(temp_df6b)[2] <- "Value"
names(temp_df7b)[2] <- "Value"
names(temp_df8b)[2] <- "Value"
names(temp_df9b)[2] <- "Value"
names(temp_df10b)[2] <- "Value"

# Agregar columna con variable categórica a cada temp_df
temp_df1a$Parameter <- "Kappa"
temp_df2a$Parameter <- "FCR"
temp_df3a$Parameter <- "NCPR"
temp_df4a$Parameter <- "Hydropathy"
temp_df5a$Parameter <- "Rg"
temp_df6a$Parameter <- "Re"
temp_df7a$Parameter <- "Aliphatic"
temp_df8a$Parameter <- "Polar"
temp_df9a$Parameter <- "Proline"
temp_df10a$Parameter <- "SCD"

# Agregar columna con variable categórica a cada temp_df
temp_df1b$Parameter <- "Kappa"
temp_df2b$Parameter <- "FCR"
temp_df3b$Parameter <- "NCPR"
temp_df4b$Parameter <- "Hydropathy"
temp_df5b$Parameter <- "Rg"
temp_df6b$Parameter <- "Re"
temp_df7b$Parameter <- "Aliphatic"
temp_df8b$Parameter <- "Polar"
temp_df9b$Parameter <- "Proline"
temp_df10b$Parameter <- "SCD"

# Unir data frames
df_sparrow1 <- rbind(temp_df1a, temp_df2a, temp_df3a, temp_df4a, temp_df5a, temp_df6a, temp_df7a, temp_df8a, temp_df9a, temp_df10a)
df_sparrow2 <- rbind(temp_df1b, temp_df2b, temp_df3b, temp_df4b, temp_df5b, temp_df6b, temp_df7b, temp_df8b, temp_df9b, temp_df10b)

# Convertir columna de parámetros en un factor
df_sparrow1$Parameter <- factor(df_sparrow1$Parameter)
df_sparrow2$Parameter <- factor(df_sparrow2$Parameter)

# Añadir columna de clasificación
df_sparrow1$Group <- "High"
df_sparrow2$Group <- "Low"

# Unir data frames df_sparrow 1 y 2
df_sparrow1 <- rbind(df_sparrow1, df_sparrow2)

# Convertir columna Group en factor
df_sparrow1$Group <- factor(df_sparrow1$Group)

# Convertir columna del grupo a factor
df_sparrow1$Group <- factor(df_sparrow1$Group)

# Anova de una vía entre los parámetros vs los grupos de alta y baja
res.aov <- aov(data = df_sparrow1, formula = Value ~ Group * Parameter)

# Verificar la estadística del anova
summary(res.aov)

# Realizar prueba post-hoc de Tukey
TukeyHSD(res.aov)

# Realizar prueba de homogeneidad de varianzas
leveneTest(data = df_sparrow1, Value ~ Group * Parameter)

# Realizar prueba de Welch (no hay homogeneidad de varianzas)
res.aov1 <- oneway.test(data = df_sparrow1, formula = Value ~ Group * Parameter)

# Encontrar diferencias entre los grupos
df_sparrow1 %>% pairwise_t_test(Value ~ Parameter)


# Constuir boxplot con los parámetros asignados
ggplot(data = df_sparrow1, aes(x = Parameter, y = Value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic() +
  labs(x = NULL, y = "CIDER value") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = "top") +
  scale_x_discrete(labels = c("Aliphatic","FCR", "Hydro", "Kappa", "NCPR", "Polar", "Proline", "Re", "Rg / 100", "SCD")) +
  scale_fill_manual(values = c("#FFA500", "#00D0D6")) +
  annotate(geom = "text", x = 1, y = 0.42, label = "ns", size = 5) +
  annotate(geom = "text", x = 2, y = 0.65, label = "***", size = 5) +
  annotate(geom = "text", x = 3, y = 0.50, label = "ns", size = 5) +
  annotate(geom = "text", x = 4, y = 0.46, label = "***", size = 5) +
  annotate(geom = "text", x = 5, y = 0.34, label = "ns", size = 5) +
  annotate(geom = "text", x = 6, y = 0.68, label = "***", size = 5) +
  annotate(geom = "text", x = 7, y = 0.25, label = "ns", size = 5) +
  annotate(geom = "text", x = 8, y = 0.57, label = "ns", size = 5) +
  annotate(geom = "text", x = 9, y = 0.54, label = "**", size = 5) +
  annotate(geom = "text", x = 10, y = 0.64, label = "**", size = 5)


# Guardar gráfico
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PDF/boxplot_parameter_cider.pdf",
       device = "pdf", width = 6, height = 4.5, units = "in", dpi = 400)

# Guardar gráfico
ggsave(filename = "E:/ME/ALL/Project_NN_proteins/PLOTS/IDR DESIGN/PNG/boxplot_parameter_cider.png",
       device = "png", width = 6, height = 4.5, units = "in", dpi = 400)


